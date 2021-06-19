The SEISCOPE optimization toolbox using Ctypes
----------------------------------------------

This repo demonstrates how it is possible to use the [SEISCOPE optimization toolbox](https://seiscope2.osug.fr/SEISCOPE-OPTIMIZATION-TOOLBOX?lang=fr) (written in Fortran) from Python. The original code is public domain and was written by Ludovic Métivier
and Romain Brossier. Minor changes to the original code have been made to allow the call of the gradient-based optimization subroutines from Python. Such changes and some improvements are listed as follows. 
 * The original source organized in 6 subdirectories, each of which is associated with one gradient-based algorithm was placed in only one folder in a modular fashion. That is, one module for each optimization algorithm grouping the procedures from each one of the old subdirectories.
 *  The Euclidean vector norm and scalar product calculations were replaced with calls to the intrinsic Fortran ```norm2``` and ```dot_product``` functions, respectively.
 *  Removing Trivial Code Duplication: same procedures in Steepest Descent and Preconditioned Nonlinear Conjugate Gradient subdirectories.
 *  Removing unused variable declarations.
 *  Vectors of lower and upper bounds (box constraints) are now optional arguments in the optimization subroutines instead of array components of a derived data type

The SEISCOPE toolbox uses a derived data type (`optim`); functionality that is not yet supported at this time by f2py - and for this reason it is used [ctypes](https://docs.python.org/3/library/ctypes.html).

The repo contains a `src` directory with the modified fortran source files and a subdirectory (`tests`). The python wrapper for the SEISCOPE optimization toolbox is found inside the `seiscope_opt_tb_wrapper` directory. Example run scripts are included in the [`examples`](https://github.com/ofmla/seiscope_opt_toolbox_w_ctypes/tree/main/seiscope_opt_tb_wrapper/examples) subdirectory within it. To build the dynamic library from the source files, use the [FoBiS](https://github.com/szaghi/FoBiS) configuration file `seiscope_opt_toolbox_w_ctypes.fobis`. 
To do this, type the following 
```
FoBiS.py build -f seiscope_opt_toolbox_w_ctypes.fobis -mode shared-gnu
```
It is also possible to build the library in static form and the `tests` programs by selecting the right `mode` flag. The full set of modes are: `static-gnu`, `static-gnu-debug`, `static-intel`, `static-intel-debug`, `shared-gnu`, `shared-gnu-debug`, `shared-intel`, `shared-intel-debug`, `tests-gnu`, `tests-gnu-debug`, `tests-intel`, `tests-intel-debug`. The modes should be self-explanatory: shared, static and tests are the modes for building (in realese, optimized form) the shared and static versions of the library and the test programs, respectively. The other 3 modes are the same, but in debug form instead of release one. -gnu use the GNU gfortran compiler while -intel the Intel one.

The examples can be run after the shared object has been built to test that wrapper works as it should. Note that you must have [Devito](https://www.devitoproject.org/) installed to run the examples. A python script `plot_curves.py` is also provide in the `python` directory. Please note that it may not be the best implementation and is intended for illustrative purposes only. The following figures were obtained with the `plot_curves.py` script after ran one of the examples.

<img src="https://github.com/ofmla/seiscope_opt_toolbox_w_ctypes/blob/main/src/python/computationalcost_curves.svg" width="425"/> <img src="https://github.com/ofmla/seiscope_opt_toolbox_w_ctypes/blob/main/src/python/convergence_curves.svg" width="425"/> 

See also
------
 * SEISCOPE optimization toolbox paper: https://library.seg.org/doi/10.1190/geo2015-0031.1
 * Original source code: http://seiscope2.osug.fr/IMG/tgz/TOOLBOX_22_09_2014_OPTIMIZATION.tgz
