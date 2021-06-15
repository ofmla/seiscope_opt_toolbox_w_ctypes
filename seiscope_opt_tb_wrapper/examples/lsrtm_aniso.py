#!/usr/bin/env python
# coding: utf-8

import pathmagic  # noqa
import numpy as np

from devito import configuration
from examples.seismic import demo_model, AcquisitionGeometry, Receiver
from examples.seismic.tti import AnisotropicWaveSolver
from devito import Function, TimeFunction, gaussian_smooth

from ctypes import c_int, c_float, c_bool

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

from interface import sotb_wrapper

'''
In this script we show a simple numerical examples for LSRTM
implemented using Devito software and the SEISCOPE optimization
toolbox (sotb) wrapper.
'''

configuration['log-level'] = 'WARNING'

# Define true and initial model
shape = (101, 101)    # Number of grid point (nx, nz)
spacing = (10., 10.)  # Grid spacing in m. The domain size is now 1km by 1km
origin = (0., 0.)     # What is the location of the top left corner.
nbl = 50
space_order = 4
dtype = np.float32

# Create solver from preset
model = demo_model('layers-tti', space_order=space_order, shape=shape, nbl=nbl,
                   dtype=dtype, spacing=spacing, fs=False, nlayers=2, vp_bottom=2.5)

model0 = demo_model('layers-tti', space_order=space_order, shape=shape, nbl=nbl,
                    dtype=dtype, spacing=spacing, fs=False, nlayers=2, vp_bottom=2.5)
gaussian_smooth(model0.vp, sigma=(5, 5))
dm = (model.vp.data**(-2) - model0.vp.data**(-2))

nshots = 11
nreceivers = 101
t0 = 0.
tn = 1000.
f0 = 0.010

# First, position source centrally in all dimensions, then set depth
src_coordinates = np.empty((1, 2))
src_coordinates[0, :] = np.array(model.domain_size) * .5
src_coordinates[0, -1] = 20.  # Depth is 20m

# Define acquisition geometry: receivers

# Initialize receivers for synthetic and imaging data
rec_coordinates = np.empty((nreceivers, 2))
rec_coordinates[:, 0] = np.linspace(0, model.domain_size[0], num=nreceivers)
rec_coordinates[:, 1] = 40.

# Geometry
geometry = AcquisitionGeometry(model, rec_coordinates, src_coordinates,
                               t0, tn, f0=f0, src_type='Ricker')

# Prepare the varying source locations
source_locations = np.empty((nshots, 2), dtype=np.float32)
source_locations[:, 0] = np.linspace(0., model.domain_size[0], num=nshots)
source_locations[:, 1] = 20.

# Compute synthetic data with forward operator
solver = AnisotropicWaveSolver(model, geometry, space_order=space_order)

# Generates true data once
dobs = np.empty((geometry.nt, nreceivers, nshots), dtype=np.float32)
usmo = np.empty((geometry.nt, model.grid.shape[0],
                 model.grid.shape[1], nshots), dtype=np.float32)
vsmo = np.empty((geometry.nt, model.grid.shape[0],
                 model.grid.shape[1], nshots), dtype=np.float32)


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

    # Create symbols to hold the gradient and residual
    grad = Function(name="grad", grid=model.grid)
    rfl = Function(name="rfl", grid=model.grid)
    rfl.data[:, :] = np.reshape(x, model.grid.shape)

    u = TimeFunction(name='u', grid=model.grid, time_order=2, space_order=space_order,
                     save=geometry.nt)
    v = TimeFunction(name='v', grid=model.grid, time_order=2, space_order=space_order,
                     save=geometry.nt)

    residual = Receiver(name='rec', grid=model.grid,
                        time_range=geometry.time_axis,
                        coordinates=geometry.rec_positions)
    objective = 0.
    grad.data[:] = 0.

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

        objective += .5*np.linalg.norm(residual.data.flatten())**2
        solver.jacobian_adjoint(rec=residual, u0=u, v0=v, vp=model0.vp, dm=grad)

    return c_float(objective), grad.data.flatten().astype(c_float)


words = ['PSTD', 'PNLCG', 'LBFGS']
a_3d_array = np.zeros((model.grid.shape[0], model.grid.shape[1], 3))

# Create an instance of the SEISCOPE optimization toolbox (sotb) Class.
sotb = sotb_wrapper()

for i, word in enumerate(words):

    # parameter initialization
    n = c_int(model.grid.shape[0]*model.grid.shape[1])  # dimension
    flag = c_int(0)                                     # first flag
    sotb.udf.conv = c_float(1e-8)   # tolerance for the stopping criterion
    sotb.udf.print_flag = c_int(1)  # print info in output files
    sotb.udf.debug = c_bool(False)  # level of details for output files
    sotb.udf.niter_max = c_int(25)  # maximum iteration number
    sotb.udf.nls_max = c_int(30)    # max number of linesearch iteration
    sotb.udf.l = c_int(5)

    # Print the derived type.
    print('Hello from Python!')
    print(sotb.udf)

    # intial guess
    X = np.zeros(model.grid.shape, dtype=c_float).reshape(-1)

    # computation of the cost and gradient associated
    # with the initial guess
    fcost, grad = lsm_gradient(X)

    # copy of grad in grad_preco: no preconditioning in
    # this test
    grad_preco = np.copy(grad)

    # optimization loop: while convergence not reached or
    # linesearch not failed, iterate

    while (flag.value != 2 and flag.value != 4):
        if word == 'PSTD':
            sotb.PSTD(n, X, fcost, grad, grad_preco, flag)
        elif word == 'PNLCG':
            sotb.PNLCG(n, X, fcost, grad, grad_preco, flag)
        else:
            sotb.LBFGS(n, X, fcost, grad, flag)
        if (flag.value == 1):
            # compute cost and gradient at point x
            fcost, grad = lsm_gradient(X)
            # no preconditioning in this test: simply copy grad in
            # grad_preco
            if word != 'LBFGS':
                grad_preco = np.copy(grad)

    # Helpful console writings
    print('END OF TEST')
    print('FINAL iterate is : ', X)
    if word == 'LBFGS':
        print('See the convergence history in iterate_'+word[:2]+'.dat')
    elif word == 'PNLCG':
        print('See the convergence history in iterate_'+word[3:]+'.dat')
    else:
        print('See the convergence history in iterate_'+word[1:-1]+'.dat')

    a_3d_array[:, :, i] = np.reshape(X, model.grid.shape)


vmax = np.amax(dm[nbl:shape[0]+nbl, nbl:shape[1]+nbl])
vmin = np.amin(dm[nbl:shape[0]+nbl, nbl:shape[1]+nbl])

mpl.rcParams['font.size'] = 8.5
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
fig.subplots_adjust(wspace=0.2)
fig.subplots_adjust(hspace=0.2)
#
im1 = ax1.imshow(dm[nbl:shape[0]+nbl, nbl:shape[1]+nbl].T, cmap=plt.cm.seismic,
                 vmin=vmin, vmax=vmax)
ax1_divider = make_axes_locatable(ax1)
cax1 = ax1_divider.append_axes("right", size="7%", pad="2%")
cb1 = plt.colorbar(im1, cax=cax1)
cb1.ax.tick_params(labelsize=8)
#
im2 = ax2.imshow(a_3d_array[nbl:shape[0]+nbl, nbl:shape[1]+nbl, 0].T,
                 cmap=plt.cm.seismic, vmin=vmin, vmax=vmax)
ax2_divider = make_axes_locatable(ax2)
cax2 = ax2_divider.append_axes("right", size="7%", pad="2%")
cb2 = plt.colorbar(im2, cax=cax2)
cb2.ax.tick_params(labelsize=8)
#
im3 = ax3.imshow(a_3d_array[nbl:shape[0]+nbl, nbl:shape[1]+nbl, 1].T,
                 cmap=plt.cm.seismic, vmin=vmin, vmax=vmax)
ax3_divider = make_axes_locatable(ax3)
cax3 = ax3_divider.append_axes("right", size="7%", pad="2%")
cb3 = plt.colorbar(im3, cax=cax3)
cb3.ax.tick_params(labelsize=8)
#
im4 = ax4.imshow(a_3d_array[nbl:shape[0]+nbl, nbl:shape[1]+nbl, 2].T,
                 cmap=plt.cm.seismic, vmin=vmin, vmax=vmax)
ax4_divider = make_axes_locatable(ax4)
cax4 = ax4_divider.append_axes("right", size="7%", pad="2%")
cb4 = plt.colorbar(im4, cax=cax4)
cb4.ax.tick_params(labelsize=8)
#
label_format = '{:,.1f}'
ticks_ylabels = (ax1.get_yticks()*0.01).tolist()
ticks_yloc = ax1.get_yticks().tolist()
ticks_xlabels = (ax1.get_xticks()*0.01).tolist()
ticks_xloc = ax1.get_xticks().tolist()
ax1.yaxis.set_major_locator(mticker.FixedLocator(ticks_yloc))
ax2.yaxis.set_major_locator(mticker.FixedLocator(ticks_yloc))
ax3.yaxis.set_major_locator(mticker.FixedLocator(ticks_yloc))
ax4.yaxis.set_major_locator(mticker.FixedLocator(ticks_yloc))

ax1.xaxis.set_major_locator(mticker.FixedLocator(ticks_xloc))
ax2.xaxis.set_major_locator(mticker.FixedLocator(ticks_xloc))
ax3.xaxis.set_major_locator(mticker.FixedLocator(ticks_xloc))
ax4.xaxis.set_major_locator(mticker.FixedLocator(ticks_xloc))

ax1.set_yticklabels([label_format.format(x) for x in ticks_ylabels])
ax2.set_yticklabels([label_format.format(x) for x in ticks_ylabels])
ax3.set_yticklabels([label_format.format(x) for x in ticks_ylabels])
ax4.set_yticklabels([label_format.format(x) for x in ticks_ylabels])

ax1.set_xticklabels([label_format.format(x) for x in ticks_xlabels])
ax2.set_xticklabels([label_format.format(x) for x in ticks_xlabels])
ax3.set_xticklabels([label_format.format(x) for x in ticks_xlabels])
ax4.set_xticklabels([label_format.format(x) for x in ticks_xlabels])

for ax in (ax1, ax2, ax3, ax4):
    ax.set(xlabel='x (km)', ylabel='Depth (km)')
for ax in (ax1, ax2, ax3, ax4):
    ax.label_outer()

plt.savefig('true_aniso.pdf')
