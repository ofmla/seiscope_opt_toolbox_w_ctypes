import ctypes
import os
import time

import numpy as np

from ctypes import POINTER, c_int, c_float, c_char_p, c_bool

# Get the location of the shared library file.
here = os.path.dirname(os.path.abspath(__file__))
lib_file = os.path.join(here, '..','fortran','lib', 'libOPTIM.so')


class UserDefined(ctypes.Structure):
    """Demonstrate how to wrap a Fortran derived type in Python using ctypes.

    Fields of the derived type are stored in the _fields_ attribute, which is a dict.
    """
    _fields_ = [
        ('debug', c_bool),
        ('threshold', c_float),
        ('print_flag', c_int),
        ('first_ls', c_bool),
        ('task', c_int),
        ('nls_max', c_int),
        ('cpt_ls', c_int),
        ('nfwd_pb', c_int),
        ('cpt_iter', c_int),
        ('niter_max', c_int),
        ('f0', c_float),
        ('fk', c_float),
        ('conv', c_float),
        ('m1', c_float),
        ('m2', c_float),
        ('mult_factor', c_float),
        ('alpha_L', c_float),
        ('alpha_R', c_float),
        ('alpha', c_float),
        ('q0', c_float),
        ('q', c_float),
        ('cpt_lbfgs',c_int),
        ('l',c_int)
    ]

    def __repr__(self):
        """Print a representation of the derived type."""
        template = (
            'UserDefined(debug={self.debug}, '
            'threshold={self.threshold}, '
            'print_flag={self.print_flag}, '
            'first_ls={self.first_ls}, '
            'task={self.task}, '
            'nls_max={self.nls_max}, '
            'cpt_ls={self.cpt_ls}, '
            'nfwd_pb={self.nfwd_pb}, '
            'cpt_iter={self.cpt_iter}, '
            'niter_max={self.niter_max}, '
            'f0={self.f0}, '
            'fk={self.fk}, '
            'conv={self.conv}, '
            'm1={self.m1}, '
            'm2={self.m2}, '
            'mult_factor={self.mult_factor}, '
            'alpha_L={self.alpha_L}, '
            'alpha_R={self.alpha_R}, '
            'alpha={self.alpha}, '
            'q0={self.q0}, '
            'conv={self.q},'
            'cpt_lbfgs={self.cpt_lbfgs}, '
            'l={self.l})'
        )
        return template.format(self=self)


def rosenbrock(X):
    """
    http://en.wikipedia.org/wiki/Rosenbrock_function
    A generalized implementation is available
    as the scipy.optimize.rosen function
    """
    a = 1. - X[0]
    b = X[1] - X[0]*X[0]
    c = X[0] - 1.
    return c_float(a*a + b*b*100.), np.array([2.*c - 400.*X[0]*b, 200.*b], dtype=c_float)


# This is how a dll/so library is loaded
lib_example = ctypes.cdll.LoadLibrary(lib_file)


def main():
    """Demonstrate how to use ctypes to use Fortran libraries."""
    # Create a UserDefined derived type in Fortran.
    conv = c_float(1e-8)
    print_flag = c_int(1)
    debug = c_bool(False)
    niter_max = c_int(10000)
    udf = UserDefined()

    # parameter initialization
    floatptr = POINTER(c_float)
    n = c_int(2)                # dimension
    flag = c_int(0)             # first flag
    udf.conv = conv             # tolerance for the stopping criterion
    udf.print_flag = print_flag # print info in output files
    udf.debug = debug           # level of details for output files
    udf.niter_max = niter_max   # maximum iteration number
    udf.nls_max = c_int(30)     # max number of linesearch iteration

    # Print the derived type.
    print('Hello from Python!')
    print(udf.__repr__())

    # intial guess
    X = np.ones(n.value, dtype=c_float)*1.50
    #grad = np.zeros(n.value, dtype=c_float)

    # computation of the cost and gradient associated
    # with the initial guess
    fcost, grad = rosenbrock(X)
    #fcost = c_float(0.)
    #lib_example.rosenbrock(X.ctypes.data_as(floatptr), ctypes.byref(fcost), grad.ctypes.data_as(floatptr))
    #print(fcost, grad)

    # copy of grad in grad_preco: no preconditioning in
    # this test
    grad_preco = np.copy(grad)

    while (flag.value != 2 and flag.value != 4):
        lib_example.PNLCG(ctypes.byref(n), X.ctypes.data_as(floatptr), ctypes.byref(fcost), grad.ctypes.data_as(floatptr),
                          grad_preco.ctypes.data_as(floatptr), ctypes.byref(udf), ctypes.byref(flag), None, None)
        #lib_example.LBFGS(ctypes.byref(n), X.ctypes.data_as(floatptr), ctypes.byref(c_float(fcost)), grad.ctypes.data_as(floatptr),
        #                  ctypes.byref(udf), ctypes.byref(flag))
        if (flag.value == 1):
            # compute cost and gradient at point x
            fcost, grad = rosenbrock(X)
            #lib_example.rosenbrock(X.ctypes.data_as(floatptr), ctypes.byref(fcost), grad.ctypes.data_as(floatptr))
            # no preconditioning in this test: simply copy grad in
            # grad_preco
            grad_preco = np.copy(grad)

    # Helpful console writings
    print('END OF TEST')
    print('FINAL iterate is : ', X)
    print('See the convergence history in iterate_CG.dat')


if __name__ == '__main__':
    main()
