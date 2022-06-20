![Build Status](https://github.com/ofmla/seiscope_opt_toolbox_w_ctypes/actions/workflows/CI.yml/badge.svg)
[![last-commit](https://img.shields.io/github/last-commit/ofmla/seiscope_opt_toolbox_w_ctypes)](https://github.com/ofmla/seiscope_opt_toolbox_w_ctypes/commits/main)
[![DOI](https://zenodo.org/badge/341695999.svg)](https://zenodo.org/badge/latestdoi/341695999)


sotb-wrapper: A Python wrapper for the SEISCOPE optimization toolbox (using Ctypes)
----------------------------------------------

This repo demonstrates how it is possible to use the [SEISCOPE optimization toolbox](https://seiscope2.osug.fr/SEISCOPE-OPTIMIZATION-TOOLBOX?lang=fr) (written in Fortran) from Python. The original code is public domain and was written by Ludovic Métivier
and Romain Brossier. Minor changes to the original code have been made to allow the call of the gradient-based optimization subroutines from Python. Such changes and some improvements are listed as follows. 
 * The original source organized in 6 subdirectories (each of which is associated with one gradient-based algorithm) was placed in only one folder in a modular fashion. That is, one module for each optimization algorithm grouping the procedures from each one of the old subdirectories.
 *  The Euclidean vector norm and scalar product calculations were replaced with calls to the intrinsic Fortran ```norm2``` and ```dot_product``` functions, respectively.
 *  Removing Trivial Code Duplication: i.e., same procedures in Steepest Descent and Preconditioned Nonlinear Conjugate Gradient subdirectories.
 *  Removing unused variable declarations.
 *  Vectors of lower and upper bounds (box constraints) are now optional arguments in the optimization subroutines instead of array components of a derived data type

The SEISCOPE toolbox uses a derived data type (`optim`); functionality that is not yet supported at this time by f2py - and for this reason it is used [ctypes](https://docs.python.org/3/library/ctypes.html). The `optim` data type is maintained, but without allocatable arrays.

The repo contains a `src` directory with the modified fortran source files and a subdirectory `tests` where each method is used to find the minimum of the banana Rosenbrock function. The python wrapper for the SEISCOPE optimization toolbox is found inside the `sotb_wrapper` directory. A `tests` directory (different from than one within `src`) includes a script to check that the wrapper has suceeded in reproducing bit-for-bit the results of the original fortran code.

Compiling
-----

To build the dynamic library from the source files, use the [FoBiS](https://github.com/szaghi/FoBiS) configuration file `seiscope_opt_toolbox_w_ctypes.fobis`. 
To do this, type the following 
```
FoBiS.py build -f seiscope_opt_toolbox_w_ctypes.fobis -mode shared-gnu
```
It is also possible to build the library in static form and the `tests` programs by changing the `mode` flag. The full set of modes are: `static-gnu`, `static-gnu-debug`, `static-intel`, `static-intel-debug`, `shared-gnu`, `shared-gnu-debug`, `shared-intel`, `shared-intel-debug`, `tests-gnu`, `tests-gnu-debug`, `tests-intel`, `tests-intel-debug`. The modes should be self-explanatory: shared, static and tests are the modes for building (in realese, optimized form) the shared and static versions of the library and the test programs, respectively. The other 3 modes are the same, but in debug form instead of release one. -gnu use the GNU gfortran compiler while -intel the Intel one.

A simple bash script `build.sh` is also provided for building the seiscope optimization toolbox with gfortran using [FoBiS](https://github.com/szaghi/FoBiS). It creates the library in static form, as originally designed as well as the tests programs. This option can be used if the user is only interested in the fortran librray.

Install sotb-wrapper
-----

After cloning the repo and build the dynamic library as indicated above in compiling step, you can install `sotb-wrapper` for development with `--editable`. What should follow is
```
git clone https://github.com/ofmla/seiscope_opt_toolbox_w_ctypes.git
cd seiscope_opt_toolbox_w_ctypes
FoBiS.py build -f seiscope_opt_toolbox_w_ctypes.fobis -mode shared-gnu
pip install -e .
```

Usage
-----

Example run scripts are included in the [`sotb_wrapper/examples`](https://github.com/ofmla/seiscope_opt_toolbox_w_ctypes/tree/main/sotb_wrapper/examples) subdirectory. The examples can be run after the shared object has been built to test that wrapper works as it should. Note that you must have [Devito](https://www.devitoproject.org/) in order to be able to run the examples. A python script `plot_curves.py` is also provide in the `examples` directory. Please note that it may not be the best implementation and is intended for illustrative purposes only. The following figures were obtained with the `plot_curves.py` script after ran one of the examples (`lsrtm_aniso.py`).

<img src="./sotb_wrapper/examples/computationalcost_curves.svg" width="425"/> <img src="./sotb_wrapper/examples/convergence_curves.svg" width="425"/> 

A tutorial in the form of a Jupyter notebook (`rosenbrock.ipynb`) is also supplied. The goal of the tutorial is show you how one can use sotb-wrapper to find a minimum for a problem, which can optionally be subject to bound constraints (also called box constraints). You can also find a simple example on calling the Fortran subroutines from a C main program in the [`c_code`](https://github.com/ofmla/seiscope_opt_toolbox_w_ctypes/tree/main/sotb_wrapper/examples/c_code) directory. The example uses the Limited-memory Broyden-Fletcher-Goldfarb-Shanno (L-BFGS) algorithm to minimize the Rosenbrock's "banana function".

License
-----

sotb-wrapper is distributed under the MIT license. See the included [`LICENSE`](https://github.com/ofmla/seiscope_opt_toolbox_w_ctypes/blob/main/LICENSE.md) file for details.

See also
------
 * SEISCOPE optimization toolbox paper: https://library.seg.org/doi/10.1190/geo2015-0031.1
 * Original source code: http://seiscope2.osug.fr/IMG/tgz/TOOLBOX_22_09_2014_OPTIMIZATION.tgz
