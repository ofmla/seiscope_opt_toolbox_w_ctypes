"""LSRTM example."""
from devito import configuration, Function, gaussian_smooth, TimeFunction
from examples.seismic import AcquisitionGeometry, demo_model, Receiver
from examples.seismic.tti import AnisotropicWaveSolver
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
import numpy as np

from sotb_wrapper import interface

"""
In this script we show a simple numerical example for LSRTM
implemented using Devito software and the SEISCOPE optimization
toolbox (sotb) wrapper.
"""

configuration["log-level"] = "WARNING"

# Define true and initial model
shape = (101, 101)  # Number of grid point (nx, nz)
spacing = (10.0, 10.0)  # Grid spacing in m. The domain size is now 1km by 1km
origin = (0.0, 0.0)  # What is the location of the top left corner.
nbl = 50
space_order = 4
dtype = np.float32

# Create solver from preset
model = demo_model(
    "layers-tti",
    space_order=space_order,
    shape=shape,
    nbl=nbl,
    dtype=dtype,
    spacing=spacing,
    fs=False,
    nlayers=2,
    vp_bottom=2.5,
)

model0 = demo_model(
    "layers-tti",
    space_order=space_order,
    shape=shape,
    nbl=nbl,
    dtype=dtype,
    spacing=spacing,
    fs=False,
    nlayers=2,
    vp_bottom=2.5,
)
gaussian_smooth(model0.vp, sigma=(5, 5))
dm = model.vp.data ** (-2) - model0.vp.data ** (-2)

nshots = 11
nreceivers = 101
t0 = 0.0
tn = 1000.0
f0 = 0.010

# First, position source centrally in all dimensions, then set depth
src_coordinates = np.empty((1, 2))
src_coordinates[0, :] = np.array(model.domain_size) * 0.5
src_coordinates[0, -1] = 20.0  # Depth is 20m

# Define acquisition geometry: receivers

# Initialize receivers for synthetic and imaging data
rec_coordinates = np.empty((nreceivers, 2))
rec_coordinates[:, 0] = np.linspace(0, model.domain_size[0], num=nreceivers)
rec_coordinates[:, 1] = 40.0

# Geometry
geometry = AcquisitionGeometry(
    model,
    rec_coordinates,
    src_coordinates,
    t0,
    tn,
    f0=f0,
    src_type="Ricker",
)

# Prepare the varying source locations
source_locations = np.empty((nshots, 2), dtype=np.float32)
source_locations[:, 0] = np.linspace(0.0, model.domain_size[0], num=nshots)
source_locations[:, 1] = 20.0

# Compute synthetic data with forward operator
solver = AnisotropicWaveSolver(model, geometry, space_order=space_order)

# Generates true data once
dobs = np.empty((geometry.nt, nreceivers, nshots), dtype=np.float32)
usmo = np.empty(
    (geometry.nt, model.grid.shape[0], model.grid.shape[1], nshots),
    dtype=np.float32,
)
vsmo = np.empty(
    (geometry.nt, model.grid.shape[0], model.grid.shape[1], nshots),
    dtype=np.float32,
)


for i in range(nshots):
    # Update source location
    geometry.src_positions[0, :] = source_locations[i, :]

    # Generate synthetic data from true model
    true_d = solver.jacobian(dm, vp=model0.vp)[0]

    # Compute smooth data and full forward wavefield u0
    u0, v0 = solver.forward(vp=model0.vp, save=True)[1:-1]

    dobs[:, :, i] = true_d.data[:, :]
    usmo[:, :, :, i] = u0.data[:, :, :]
    vsmo[:, :, :, i] = v0.data[:, :, :]


def lsm_gradient(x):
    """Computes the gradient for all shots.

    Args:
        x: Vector of Reflectivity

    Returns:
        objective: Cost function value associated with x
        grad: Gradient at x
    """
    # Create symbols to hold the gradient and residual
    grad = Function(name="grad", grid=model.grid)
    rfl = Function(name="rfl", grid=model.grid)
    rfl.data[:, :] = np.reshape(x, model.grid.shape)

    u = TimeFunction(
        name="u",
        grid=model.grid,
        time_order=2,
        space_order=space_order,
        save=geometry.nt,
    )
    v = TimeFunction(
        name="v",
        grid=model.grid,
        time_order=2,
        space_order=space_order,
        save=geometry.nt,
    )

    residual = Receiver(
        name="rec",
        grid=model.grid,
        time_range=geometry.time_axis,
        coordinates=geometry.rec_positions,
    )
    objective = 0.0
    grad.data[:] = 0.0

    for i in range(nshots):
        # Update source location
        geometry.src_positions[0, :] = source_locations[i, :]

        # Compute smooth data and full forward wavefield u0
        u.data[:, :, :] = usmo[:, :, :, i]
        v.data[:, :, :] = vsmo[:, :, :, i]

        # Compute smooth data and full forward wavefield u0
        smooth_d = solver.jacobian(rfl, vp=model0.vp)[0]

        # Compute gradient from data residual and update objective function
        residual.data[:] = smooth_d.data[:] - dobs[:, :, i]

        objective += 0.5 * np.linalg.norm(residual.data.flatten()) ** 2
        solver.jacobian_adjoint(rec=residual, u0=u, v0=v, vp=model0.vp, dm=grad)

    return objective, grad.data.flatten().astype(np.float32)


methods = ["PSTD", "PNLCG", "LBFGS"]
string = "See the convergence history in iterate_"
comments = {
    "PSTD": string + methods[0][1:-1] + ".dat",
    "PNLCG": string + methods[1][3:] + ".dat",
    "LBFGS": string + methods[2][:2] + ".dat",
}
a_3d_array = np.zeros((model.grid.shape[0], model.grid.shape[1], 3))

# Create an instance of the SEISCOPE optimization toolbox (sotb) Class.
sotb = interface.sotb_wrapper()

for i, method in enumerate(methods):

    n = model.grid.shape[0] * model.grid.shape[1]  # dimension
    flag = 0  # first flag
    # define some entries for the udf dictionary
    conv = 1e-8  # tolerance for the stopping criterion
    print_flag = 1  # print info in output files
    debug = False  # level of details for output files
    niter_max = 25  # maximum iteration number
    l = 5  # maximum number of stored pairs used for the l-BFGS approximation

    # intial guess
    X = np.zeros(model.grid.shape, dtype=np.float32).reshape(-1)

    # computation of the cost and gradient associated
    # with the initial guess
    fcost, grad = lsm_gradient(X)

    # parameter initialization
    # reset udf dict as only a few itens are initialized by `set_inputs` function
    if i > 0:
        sotb.udf = interface.UserDefined()
    sotb.set_inputs(
        fcost, niter_max, conv=conv, print_flag=print_flag, l=l, debug=debug
    )

    # Print the derived type.
    print(sotb.udf)

    # copy of grad in grad_preco: no preconditioning in
    # this test
    if method != "LBFGS":
        grad_preco = np.copy(grad)

    # optimization loop: while convergence not reached or
    # linesearch not failed, iterate

    while flag != 2 and flag != 4:
        if method == "PSTD":
            flag = sotb.PSTD(n, X, fcost, grad, grad_preco, flag)
        elif method == "PNLCG":
            flag = sotb.PNLCG(n, X, fcost, grad, grad_preco, flag)
        else:
            flag = sotb.LBFGS(n, X, fcost, grad, flag)
        if flag == 1:
            # compute cost and gradient at point x
            fcost, grad = lsm_gradient(X)
            # no preconditioning in this test: simply copy grad in
            # grad_preco
            if method != "LBFGS":
                grad_preco = np.copy(grad)

    # Helpful console writings
    print("END OF TEST")
    print("FINAL iterate is : ", X)
    print(comments[method])

    a_3d_array[:, :, i] = np.reshape(X, model.grid.shape)

# Plot results
vmax = np.amax(dm[nbl : shape[0] + nbl, nbl : shape[1] + nbl])
vmin = np.amin(dm[nbl : shape[0] + nbl, nbl : shape[1] + nbl])

mpl.rcParams["font.size"] = 8.5
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
fig.subplots_adjust(wspace=0.2)
fig.subplots_adjust(hspace=0.2)
label_format = "{:,.1f}"
#
for count, ax in enumerate([ax1, ax2, ax3, ax4]):
    if count == 0:
        im = ax.imshow(
            dm[nbl : shape[0] + nbl, nbl : shape[1] + nbl].T,
            cmap=plt.cm.seismic,
            vmin=vmin,
            vmax=vmax,
        )
    else:
        im = ax.imshow(
            a_3d_array[nbl : shape[0] + nbl, nbl : shape[1] + nbl, count - 1].T,
            cmap=plt.cm.seismic,
            vmin=vmin,
            vmax=vmax,
        )
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

plt.savefig("true_aniso.pdf")
