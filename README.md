The SEISCOPE optimization toolbox using Ctypes
----------------------------------------------

This repo demonstrates how it is possible to use the [SEISCOPE optimization toolbox](https://seiscope2.osug.fr/SEISCOPE-OPTIMIZATION-TOOLBOX?lang=fr) (written in Fortran) from Python. The original code is public domain and was written by Ludovic Métivier
and Romain Brossier. Minor changes to the original code have been made to allow the call of the gradient-based optimization subroutines from Python. Such changes are and some improvements are listed as follows. 
 * The original source organized in 6 subdirectories, each of which is associated with one gradient-based algorithm was placed in only one folder in a modular fashion. That is, one module for each optimization algorithm grouping the procedures from each one of the old subdirectories.
 *  The Euclidean vector norm and scalar product calculations were replaced with calls to the intrinsic Fortran ```norm2``` and ```dot_product``` functions, respectively.
 *  New Makefile for building a dynamic version of the library with a debug option.
 *  Removing Trivial Code Duplication: same procedures in Steepest Descent and Preconditioned Nonlinear Conjugate Gradient subdirectories.
 *  Removing unused variable declarations.

The SEISCOPE toolbox uses a derived data type (`optim`); functionality that is not yet supported at this time by f2py - and for this reason it is used [ctypes](https://docs.python.org/3/library/ctypes.html).

The repo contains a single `src` directory with two subdirectories `fortran` and `python`. The modified source files along the `Makefile` are in the `fortran ` folder, while `python` contains two python scripts with examples where the toolbox is used: The *Rosenbrock function* and *Least Square Reverse Time Migration (LSQRTM)* optimization problems. To build the dynamic library from the source files, execute the `Makefile` in the `fortran` directory. The examples can be run after the shared
object has been built to test that it works as it should. Note that you must have [Devito](https://www.devitoproject.org/) installed to run the *LSQRTM* example.

A python script `plot_curves.py` is also provide in the `python` directory. Please note that it may not be the best implementation and is intended for illustrative purposes only. The following figures were obtained after the LSQRTM run with the `plot_curves.py` script.

<img src="https://github.com/ofmla/seiscope_opt_toolbox_w_ctypes/blob/main/src/python/computationalcost_curves.svg" width="425"/> <img src="https://github.com/ofmla/seiscope_opt_toolbox_w_ctypes/blob/main/src/python/convergence_curves.svg" width="425"/> 

See also
------
 * SEISCOPE optimization toolbox paper: https://library.seg.org/doi/10.1190/geo2015-0031.1
 * Original source code: http://seiscope2.osug.fr/IMG/tgz/TOOLBOX_22_09_2014_OPTIMIZATION.tgz
