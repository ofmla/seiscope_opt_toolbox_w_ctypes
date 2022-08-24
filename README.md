[![GitHub release (latest by date)](https://img.shields.io/github/v/release/ofmla/seiscope_opt_toolbox_w_ctypes)](https://github.com/ofmla/seiscope_opt_toolbox_w_ctypes/releases/tag/v1.0.1)
![Build Status](https://github.com/ofmla/seiscope_opt_toolbox_w_ctypes/actions/workflows/python.yml/badge.svg)
[![last-commit](https://img.shields.io/github/last-commit/ofmla/seiscope_opt_toolbox_w_ctypes)](https://github.com/ofmla/seiscope_opt_toolbox_w_ctypes/commits/main)
[![DOI](https://zenodo.org/badge/341695999.svg)](https://zenodo.org/badge/latestdoi/341695999)
[![codecov](https://codecov.io/gh/ofmla/seiscope_opt_toolbox_w_ctypes/branch/main/graph/badge.svg?token=BN7NK8A9OJ)](https://codecov.io/gh/ofmla/seiscope_opt_toolbox_w_ctypes)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)


# sotb-wrapper: A Python wrapper for the SEISCOPE optimization toolbox (using Ctypes)

This repo demonstrates how it is possible to use the [SEISCOPE optimization toolbox](https://seiscope2.osug.fr/SEISCOPE-OPTIMIZATION-TOOLBOX?lang=fr) (written in Fortran) from Python. The original code is public domain and was written by Ludovic Métivier
and Romain Brossier. Minor changes to the original code have been made to allow the call of the gradient-based optimization subroutines from Python. Such changes and some improvements are listed as follows. 
 * The original source organized in 6 subdirectories (each of which is associated with one gradient-based algorithm) was placed in only one folder in a modular fashion. That is, one module for each optimization algorithm grouping the procedures from each one of the old subdirectories.
 *  The Euclidean vector norm and scalar product calculations were replaced with calls to the intrinsic Fortran ```norm2``` and ```dot_product``` functions, respectively.
 *  Removing Trivial Code Duplication: i.e., same procedures in Steepest Descent and Preconditioned Nonlinear Conjugate Gradient subdirectories.
 *  Removing unused variable declarations.
 *  Vectors of lower and upper bounds (box constraints) are now optional arguments in the optimization subroutines instead of array components of a derived data type

The SEISCOPE toolbox uses a derived data type (`optim`); functionality that is not yet supported at this time by f2py - and for this reason it is used [ctypes](https://docs.python.org/3/library/ctypes.html). The `optim` data type is maintained, but without allocatable arrays.

The repo contains a `src` directory with the modified fortran source files and another named `apps` where each method is used to find the minimum of the banana Rosenbrock function. The python wrapper for the SEISCOPE optimization toolbox is found inside the `sotb_wrapper` directory. A `test` directory includes a script to check that the wrapper has suceeded in reproducing the results of the original fortran code.

## Install Seiscope optimization toolbox (sotb)

If you only want to use the Fortran library, you can simply clone the repo and build it with [Fortran Package Manager](https://github.com/fortran-lang/fpm) or [CMake](https://cmake.org/). In the first case you just need to run 
```bash
fpm build --profile release
```
This command creates the library in static form, as originally designed as well as executable files from demo codes. 

To use `sotb` within your fpm project, add the following to your `fpm.toml` file:

```yml
[dependencies]
sotb = { git="https://github.com/ofmla/seiscope_opt_toolbox_w_ctypes.git" }
```
In the second case, you can run a workflow as the following:
```bash
FC=gfortran cmake -B _build -DCMAKE_INSTALL_PREFIX=$PREFIX -DCMAKE_BUILD_TYPE=Release
cmake --build _build
cmake --install _build
```
where you need to replace `$PREFIX` with the desired directory.

Examples of use of `sotb` can be found in the `app` folder, which contains a folder with an example for each one of the optimization algorithms available in the library. The executable files for each example are built with `cmake` invocation above and made available at `$PREFIX/bin` folder. As mentioned before, when you use `fpm`, executable files for the examples are also created. In this latest case, you can use `fpm run --profile release <test_name>` to run an specific example. So, if you want to run the example that uses the limited-memory version of Broyden-Fletcher-Goldfarb-Shanno (L-BFGS) algorithm, simply run `fpm run --profile release test_LBFGS`. If you run `fpm run --profile release` you can see the names of the six available examples. You can also find a simple example on calling the Fortran subroutines from a C main program in the [`c_code`](https://github.com/ofmla/seiscope_opt_toolbox_w_ctypes/tree/main/examples/c_code) directory. The example uses the L-BFGS to minimize the Rosenbrock's "banana function". Assuming that `$PREFIX` points to the repository root directory, you can create the executable from [`c_code`](https://github.com/ofmla/seiscope_opt_toolbox_w_ctypes/tree/main/examples/c_code) directory, by running 
```bash
cmake -S. -B _build -DCMAKE_PREFIX_PATH="`pwd`/../../"
cmake --build _build
```
or
```bash
cmake -S. -B _build -Dsotb_DIR="`pwd`/../../lib/cmake/sotb"
cmake --build _build
```
## Install the python wrapper of sotb (sotb-wrapper)

To install the Python API with the embedded `sotb` shared library you can use pip.
```bash
pip install sotb-wrapper
```
It is also possible to install directly from the GitHub repository. You will need a Fortran compiler such as GFortran to compile the shared library but the whole process is automated via [scikit-build](https://github.com/scikit-build/scikit-build).
```bash
pip install git+https://github.com/ofmla/seiscope_opt_toolbox_w_ctypes
```
### Usage

The following example demonstrates how to define and solve the classical Rosenbrock problem
```python
import numpy as np
from sotb_wrapper import interface

# Declare the objective function
def rosenbrock(X):
    """
    http://en.wikipedia.org/wiki/Rosenbrock_function
    A generalized implementation is available
    as the scipy.optimize.rosen function
    """
    a = 1. - X[0]
    b = X[1] - X[0]*X[0]
    return a*a + b*b*100., np.array([-a*2. - 400.*X[0]*b, 200.*b], dtype=np.float32)
    
# Create an instance of the SEISCOPE optimization toolbox wrapper (sotb_wrapper) Class. 
sotb = interface.sotb_wrapper()
n = 2 # dimension
flag = 0 # first flag; 0 means initialization
X = np.ones(2, dtype=np.float32)*-1. # initial guess

# computation of the cost and gradient associated with the initial guess
fcost, grad = rosenbrock(X)
# copy of grad in grad_preco: no preconditioning in this test
grad_preco = np.copy(grad)
# Set parameters of the UserDefined derived type in Fortran (ctype structure).
# The first two parameters are mandatory; all others are optional. 
sotb.set_inputs(fcost, 10000, conv=1e-8, l=10)

# optimization loop: while convergence not reached or linesearch not failed, iterate
while (flag != 2 and flag != 4):
    flag = sotb.PSTD(n, X, fcost, grad, grad_preco, flag)
    if flag == 1:
        # compute cost and gradient at point x
        fcost, grad = rosenbrock(X)
        # no preconditioning in this test: simply copy grad in grad_preco
        grad_preco = np.copy(grad)
print('FINAL iterate is : ', X)
```

The code above is part of a tutorial in the form of a Jupyter notebook (`rosenbrock.ipynb`) provided in the [`examples`](https://github.com/ofmla/seiscope_opt_toolbox_w_ctypes/tree/main/examples) subdirectory. The goal of the tutorial is show you how one can use sotb-wrapper to find a minimum for a problem, which can optionally be subject to bound constraints (also called box constraints). The directory also includes examples in the context of geophysical inversion. Note that you must have [Devito](https://www.devitoproject.org/) in order to be able to run them. A python script `plot_curves.py` is also provide in the `examples` directory. It may not be the best implementation and is intended for illustrative purposes only. 

The following figures were obtained with the `plot_curves.py` script after ran one of the examples (`lsrtm_aniso.py`).

<img src="./examples/computationalcost_curves.svg" width="425"/> <img src="./examples/convergence_curves.svg" width="425"/> 

License
-----

sotb-wrapper is distributed under the MIT license. See the included [`LICENSE`](https://github.com/ofmla/seiscope_opt_toolbox_w_ctypes/blob/main/LICENSE.md) file for details.

See also
------
 * SEISCOPE optimization toolbox paper: https://library.seg.org/doi/10.1190/geo2015-0031.1
 * Original source code: http://seiscope2.osug.fr/IMG/tgz/TOOLBOX_22_09_2014_OPTIMIZATION.tgz
