The SEISCOPE optimization toolbox using Ctypes
----------------------------------------------

This repo demonstrates how it is possible to use the SEISCOPE optimization toolbox (written in Fortran) from Python. Minor changes to the original code have been made to allow the call of the gradient-based optimization subroutines from Python. The SEISCOPE toolbox uses a derived data type (`optim`); functionality that is not yet supported at this time by f2py - and for this reason it is used ctypes.

The repo contains a single `src` directory with two subdirectories `fortran` and `python`. The modified source files along the `Makefile` are in the `fortran ` folder, while `python` contains two examples where the toolbox is used: The *Rosenbrock function* and *Least Square Reverse Time Migration (LSQRTM)* optimization problems. To build the dynamic library from the source files, execute the `Makefile` in the `fortran` directory. The examples can be run after the shared
object has been built to test that it works as it should. Note that you must have `Devito` installed to run the *LSQRTM* example.

A python script `plot_curves.py` is also provide in the `python` directory. The following figures were obtained after the LSQRTM run, with the `plot_curves.py` script.

<img src="https://github.com/ofmla/seiscope_opt_toolbox_w_ctypes/blob/main/src/python/computationalcost_curves.svg" width="460"/> <img src="https://github.com/ofmla/seiscope_opt_toolbox_w_ctypes/blob/main/src/python/convergence_curves.svg" width="460"/> 
