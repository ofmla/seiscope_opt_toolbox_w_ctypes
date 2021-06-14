import numpy as np
import ctypes
from interface import sotb_wrapper, lib_sotb
from ctypes import c_int, c_float, c_bool, POINTER
from os import path


def rosenbrock(X):
    """
    http://en.wikipedia.org/wiki/Rosenbrock_function
    A generalized implementation is available
    as the scipy.optimize.rosen function
    """
    a = 1. - X[0]
    b = X[1] - X[0]**2
    c = X[0] - 1.
    return c_float(a**2 + 100.*b**2), np.array([2.*c - 400.*X[0]*b, 200.*b], dtype=np.float32)


def rosenbrock_hess(X, d):

    a = (1200.*X[0]**2 - 400.*X[1] + 2.)*d[0] - 400*X[0]*d[1]
    b = -400.*X[0]*d[0] + 200.*d[1]
    return np.array([a, b], dtype=np.float32)


def main():
    """
    Demonstrate how to use the SEISCOPE optimization toolbox (sotb) wrapper to
    find the minimum of Rosenbrock's banana function (a famous test case for
    optimization software)
    """
    # Create an instance of the SEISCOPE optimization toolbox (sotb) Class.
    floatptr = POINTER(c_float)
    sotb = sotb_wrapper()

    methods = ['PSTD', 'PNLCG', 'LBFGS', 'TRN']

    for i, method in enumerate(methods):
        # Set some fields of the UserDefined derived type in Fortran (ctype structure).
        # parameter initialization
        n = c_int(2)                       # dimension
        flag = c_int(0)                    # first flag
        sotb.udf.conv = c_float(1e-8)      # tolerance for the stopping criterion
        sotb.udf.print_flag = c_int(1)     # print info in output files
        sotb.udf.debug = c_bool(False)     # level of details for output files
        if method == 'TRN':
            sotb.udf.niter_max = c_int(100)  # maximum iteration number
        else:
            sotb.udf.niter_max = c_int(10000)
        sotb.udf.nls_max = c_int(30)       # max number of linesearch iteration
        sotb.udf.l = c_int(20)
        sotb.udf.niter_max_CG = c_int(5)   # max no. of inner CG iterations

        # intial guess
        X = np.ones(2, dtype=np.float32)*1.5

        # computation of the cost and gradient associated
        # with the initial guess
        # fcost, grad = rosenbrock(X)
        grad = np.zeros(n.value, dtype=np.float32)
        d = np.zeros(n.value, dtype=np.float32)
        Hd = np.zeros(n.value, dtype=np.float32)
        fcost = c_float(0.)
        lib_sotb.rosenbrock(X.ctypes.data_as(floatptr), ctypes.byref(fcost),
                            grad.ctypes.data_as(floatptr))

        # copy of grad in grad_preco: no preconditioning in
        # this test
        grad_preco = np.zeros(2, dtype=np.float32)
        if method == 'PNLCG':
            grad_preco = np.copy(grad)

        # optimization loop: while convergence not reached or
        # linesearch not failed, iterate

        while (flag.value != 2 and flag.value != 4):
            if method == 'PSTD':
                sotb.PSTD(n, X, fcost, grad, grad_preco, flag)
            elif method == 'PNLCG':
                sotb.PNLCG(n, X, fcost, grad, grad_preco, flag)
            elif method == 'LBFGS':
                sotb.LBFGS(n, X, fcost, grad, flag)
            else:
                sotb.TRN(n, X, fcost, grad, d, Hd, flag)
            if (flag.value == 1):
                # compute cost and gradient at point x
                # fcost, grad = rosenbrock(X)
                lib_sotb.rosenbrock(X.ctypes.data_as(floatptr), ctypes.byref(fcost),
                                    grad.ctypes.data_as(floatptr))
                # no preconditioning in this test: simply copy grad in
                # grad_preco
                if method != 'LBFGS':
                    grad_preco = np.copy(grad)
            elif (flag.value == 7):
                # compute d by the Hessian operator and store in Hd
                lib_sotb.rosenbrock_hess(X.ctypes.data_as(floatptr),
                                         d.ctypes.data_as(floatptr),
                                         Hd.ctypes.data_as(floatptr))

        # Helpful console writings
        print('END OF TEST')
        print('FINAL iterate is : ', X)
        if method == 'LBFGS':
            print('See the convergence history in iterate_'+method[:2]+'.dat')
        elif method == 'PNLCG':
            print('See the convergence history in iterate_'+method[3:]+'.dat')
        elif method == 'LBFGS':
            print('See the convergence history in iterate_'+method[1:-1]+'.dat')
        else:
            print('See the convergence history in iterate_'+method+'.dat and iterate_'+method+'_CG.dat')

    assert path.isfile('iterate_ST.dat')
    assert path.isfile('iterate_CG.dat')
    assert path.isfile('iterate_LB.dat')
    assert path.isfile('iterate_TRN.dat')

if __name__ == '__main__':
    main()

