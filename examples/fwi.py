"""FWI example."""
from devito import Function
from distributed import Client, LocalCluster, wait
from examples.seismic import AcquisitionGeometry, demo_model, Receiver
from examples.seismic.acoustic import AcousticWaveSolver
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import numpy as np

from sotb_wrapper import interface


def forward_modeling_single_shot(model, geometry, save=False, dt=4.0):
    """Serial modeling function.

    Args:
        model: instance of SeismicModel class that has physical parameter
            attributes, such as `vp`
        geometry: Acquisition geometry for a single shot
        save: Whether or not to save the entire wavefield. Defaults to False
        dt: Sampling interval used for resampling the modeled shot. Defaults to 4.0 ms

    Returns:
        d_obs: Resampled shot at intervals of `dt` ms
        u0: Entire wavefield or last three snapshots, depending on `save` value
    """
    solver = AcousticWaveSolver(model, geometry, space_order=4)
    d_obs, u0 = solver.forward(vp=model.vp, save=save)[0:2]
    return d_obs.resample(dt), u0


def forward_modeling_multi_shots(client, model, geometry, save=False, dt=4.0):
    """Parallel modeling function.

    Args:
        client: Dask scheduler
        model: instance of SeismicModel class that has physical parameter
            attributes, such as `vp`
        geometry: Acquisition geometry
        save: Whether or not to save the entire wavefield. Defaults to False
        dt: Sampling interval used for resampling the modeled shot. Defaults to 4.0 ms

    Returns:
        shots: list of modeled shots
    """
    futures = []
    for i in range(geometry.nsrc):

        # Geometry for current shot
        geometry_i = AcquisitionGeometry(
            model,
            geometry.rec_positions,
            geometry.src_positions[i, :],
            geometry.t0,
            geometry.tn,
            f0=geometry.f0,
            src_type=geometry.src_type,
        )

        # Call serial modeling function for each index
        futures.append(
            client.submit(
                forward_modeling_single_shot, model, geometry_i, save=save, dt=dt
            )
        )

    # Wait for all workers to finish and collect shots
    wait(futures)
    shots = []
    for i in range(geometry.nsrc):
        shots.append(futures[i].result()[0])

    return shots


def fwi_objective_single_shot(model, geometry, d_obs):
    """Serial FWI objective function.

    Args:
        model: instance of SeismicModel class that has physical parameter
            attributes, such as `vp`
        geometry: Acquisition geometry for a single shot
        d_obs: True data (one single shot)

    Returns:
        fval: Cost function value for one shot
        grad: Gradient for one source location, without absorbing boundaries
    """
    # Devito objects for gradient and data residual
    grad = Function(name="grad", grid=model.grid)
    residual = Receiver(
        name="rec",
        grid=model.grid,
        time_range=geometry.time_axis,
        coordinates=geometry.rec_positions,
    )
    solver = AcousticWaveSolver(model, geometry, space_order=4)

    # Predicted data and residual
    d_pred, u0 = solver.forward(vp=model.vp, save=True)[0:2]
    dobs_resampled = d_obs.resample(geometry.dt).data[:][0 : d_pred.data.shape[0], :]
    residual.data[:] = d_pred.data[:] - dobs_resampled

    # Function value and gradient
    fval = 0.5 * np.linalg.norm(residual.data.flatten()) ** 2
    solver.gradient(rec=residual, u=u0, vp=model.vp, grad=grad)

    # Convert to numpy array and remove absorbing boundaries
    grad_crop = np.array(grad.data[:])[model.nbl : -model.nbl, model.nbl : -model.nbl]

    return fval, grad_crop


def fwi_objective_multi_shots(client, model, geometry, d_obs):
    """Parallel FWI objective function.

    Args:
        client: Dask scheduler
        model: instance of SeismicModel class that has physical parameter
            attributes, such as `vp`
        geometry: Acquisition geometry
        d_obs: True data (all shots)

    Returns:
        fval: Cost function value associated with the current model
        grad: Gradient at x
    """
    futures = []
    for i in range(geometry.nsrc):

        # Geometry for current shot
        geometry_i = AcquisitionGeometry(
            model,
            geometry.rec_positions,
            geometry.src_positions[i, :],
            geometry.t0,
            geometry.tn,
            f0=geometry.f0,
            src_type=geometry.src_type,
        )

        # Call serial FWI objective function for each shot location
        futures.append(
            client.submit(fwi_objective_single_shot, model, geometry_i, d_obs[i])
        )

    # Wait for all workers to finish and collect function values and gradients
    wait(futures)
    fval = 0.0
    grad = np.zeros(model.shape)
    for i in range(geometry.nsrc):
        fval += futures[i].result()[0]
        grad += futures[i].result()[1]

    return fval, grad


def loss(c, x, model, geometry, d_obs):
    """Calculates the gradient and objective function value.

    Args:
        c: Dask scheduler.
        x: squared slowness [s^2/km^2]
        model: instance of SeismicModel class that has physical parameter
            attributes, such as `vp`
        geometry: Acquisition geometry
        d_obs: True data

    Returns:
        fval: Cost function value associated with the current model
        grad: Gradient at x
    """
    # Convert x to velocity
    v_curr = 1.0 / np.sqrt(x.reshape(model.shape))

    # Overwrite current velocity in geometry (don't update boundary region)
    model.update("vp", v_curr.reshape(model.shape))

    # Evaluate objective function
    fval, grad = fwi_objective_multi_shots(c, model, geometry, d_obs)
    return fval, grad.flatten().astype(np.float32)


def main(c):
    """A basic FWI implementation. It uses the l-BFGS method for the optimization.

    Args:
        c: Dask scheduler
    """
    # Set up velocity model
    shape = (101, 101)  # Number of grid points (nx, nz).
    spacing = (10.0, 10.0)  # Grid spacing in m. The domain size is now 1km by 1km.
    origin = (0, 0)  # Need origin to define relative source and receiver locations.
    nbl = 40

    # True model
    model1 = demo_model(
        "circle-isotropic",
        vp_circle=3.0,
        vp_background=2.5,
        origin=origin,
        shape=shape,
        spacing=spacing,
        nbl=nbl,
    )

    # Initial model
    model0 = demo_model(
        "circle-isotropic",
        vp_circle=2.5,
        vp_background=2.5,
        origin=origin,
        shape=shape,
        spacing=spacing,
        nbl=nbl,
        grid=model1.grid,
    )

    # Set up acquisiton geometry
    t0 = 0.0
    tn = 1000.0
    f0 = 0.010

    # Set up source geometry, but define 5 sources instead of just one.
    nsources = 5
    src_coordinates = np.empty((nsources, 2))
    src_coordinates[:, 1] = np.linspace(0, model1.domain_size[0], num=nsources)
    src_coordinates[:, 0] = 20.0  # Source depth is 20m

    # Initialize receivers for synthetic and imaging data
    nreceivers = 101
    rec_coordinates = np.empty((nreceivers, 2))
    rec_coordinates[:, 1] = np.linspace(
        spacing[0], model1.domain_size[0] - spacing[0], num=nreceivers
    )
    rec_coordinates[:, 0] = 980.0  # Receiver depth
    # Set up geometry objects for observed and predicted data
    geometry1 = AcquisitionGeometry(
        model1, rec_coordinates, src_coordinates, t0, tn, f0=f0, src_type="Ricker"
    )
    geometry0 = AcquisitionGeometry(
        model0, rec_coordinates, src_coordinates, t0, tn, f0=f0, src_type="Ricker"
    )

    # Compute observed data in parallel (inverse crime). In real life we would
    # read the SEG-Y data here.
    d_obs = forward_modeling_multi_shots(c, model1, geometry1, save=False)

    # Box contraints
    vmin = 1.4  # do not allow velocities slower than water
    vmax = 4.0
    lb = [1.0 / vmax**2 for _ in range(np.prod(model0.shape))]  # in [s^2/km^2]
    ub = [1.0 / vmin**2 for _ in range(np.prod(model0.shape))]  # in [s^2/km^2]
    lb = np.array(lb).astype(np.float32)
    ub = np.array(ub).astype(np.float32)

    # Create an instance of the SEISCOPE optimization toolbox (sotb) Class.
    sotb = interface.sotb_wrapper()

    n = shape[0] * shape[1]  # dimension
    flag = 0  # first flag

    # Define some fields of the UserDefined derived type in Fortran (ctype structure).
    conv = 1e-8  # tolerance for the stopping criterion
    print_flag = 1  # print info in output files
    debug = False  # level of details for output files
    niter_max = 5  # maximum iteration number
    l = 5  # maximum number of stored pairs used for the l-BFGS approximation

    # Initial guess
    v0 = model0.vp.data[model0.nbl : -model0.nbl, model0.nbl : -model0.nbl]
    X = 1.0 / (v0.reshape(-1).astype(np.float32)) ** 2

    # computation of the cost and gradient associated
    # with the initial guess
    fcost, grad = loss(c, X, model0, geometry0, d_obs)

    # parameter initialization
    sotb.set_inputs(
        fcost, niter_max, conv=conv, print_flag=print_flag, l=l, debug=debug
    )

    # Print the derived type.
    print(sotb.udf)

    while flag != 2 and flag != 4:
        flag = sotb.LBFGS(n, X, fcost, grad, flag, lb, ub)
        if flag == 1:
            # compute cost and gradient at point x
            fcost, grad = loss(c, X, model0, geometry0, d_obs)

    # Helpful console writings
    print("END OF TEST")
    print("FINAL iterate is : ", X)
    print("See the convergence history in iterate_LB.dat")

    # Plot FWI result
    vp = 1.0 / np.sqrt(X.reshape(model1.shape))
    vmin = 2.4
    vmax = 2.8

    mpl.rcParams["font.size"] = 8.5
    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.subplots_adjust(wspace=0.25)
    label_format = "{:,.1f}"
    #
    for count, ax in enumerate([ax1, ax2]):
        if count == 0:
            im = ax.imshow(
                model1.vp.data[nbl:-nbl, nbl:-nbl].T,
                cmap=plt.cm.cividis,
                vmin=vmin,
                vmax=vmax,
            )
        else:
            im = ax.imshow(vp.T, cmap=plt.cm.cividis, vmin=vmin, vmax=vmax)
        ax_divider = make_axes_locatable(ax)
        cax = ax_divider.append_axes("right", size="7%", pad="2%")
        cb = plt.colorbar(im, cax=cax)
        cb.ax.tick_params(labelsize=8)
        #
        ticks_ylabels = (ax.get_yticks() * 0.01).tolist()
        ticks_yloc = ax.get_yticks().tolist()
        ticks_xlabels = (ax.get_xticks() * 0.01).tolist()
        ticks_xloc = ax.get_xticks().tolist()
        ax.yaxis.set_major_locator(mticker.FixedLocator(ticks_yloc))
        ax.xaxis.set_major_locator(mticker.FixedLocator(ticks_xloc))
        ax.set_yticklabels([label_format.format(y) for y in ticks_ylabels])
        ax.set_xticklabels([label_format.format(x) for x in ticks_xlabels])
        ax.set(xlabel="x (km)", ylabel="Depth (km)")
        ax.label_outer()

    plt.savefig("circle_isotropic_inversion.pdf")


if __name__ == "__main__":
    r"""
    This script demonstrates how we can set up a basic FWI framework
    with gradient-based optimization algorithms from the SEISCOPE
    optimization toolbox (sotb) wrapper. The script is basically a copy
    of the devito FWI tutorial (https://github.com/devitocodes/devito/
    blob/master/examples/seismic/tutorials/04_dask.ipynb) with the addition
    of the sotb wrapper. It uses a simple toy example for validation of the
    code.
    """
    print("start")
    # Start Dask cluster
    cluster = LocalCluster(n_workers=5, death_timeout=600, asynchronous=False)
    c = Client(cluster)
    main(c)
