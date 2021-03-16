import numpy as np
from interface import sotb_wrapper
from ctypes import c_int, c_float, c_bool

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

def main():
    """
    Demonstrate how to use the SEISCOPE optimization toolbox (sotb) wrapper to
    find the minimum of Rosenbrock's banana function (a famous test case for
    optimization software)
    """
    # Create ctype variables to build the UserDefined derived type in Fortran.
    conv = c_float(1e-8)
    print_flag = c_int(1)
    debug = c_bool(False)
    niter_max = c_int(10000)
    
    # Create an instance of the SEISCOPE optimization toolbox (sotb) Class.
    sotb = sotb_wrapper()
    
	# Set some fields of the UserDefined derived type in Fortran (ctype structure).
    # parameter initialization
    n = c_int(2)                     # dimension
    flag = c_int(0)                  # first flag
    sotb.udf.conv = conv             # tolerance for the stopping criterion
    sotb.udf.print_flag = print_flag # print info in output files
    sotb.udf.debug = debug           # level of details for output files
    sotb.udf.niter_max = niter_max   # maximum iteration number
    sotb.udf.nls_max = c_int(30)     # max number of linesearch iteration

    # Print the derived type.
    print('Hello from Python!')
    print(sotb.udf)

    # intial guess
    X = np.ones(n.value, dtype=c_float)*1.50

    # computation of the cost and gradient associated
    # with the initial guess
    fcost, grad = rosenbrock(X)

    # copy of grad in grad_preco: no preconditioning in
    # this test
    grad_preco = np.copy(grad)

    while (flag.value != 2 and flag.value != 4):
        sotb.PNLCG(n, X, fcost, grad, grad_preco, flag)
        if (flag.value == 1):
            # compute cost and gradient at point x
            fcost, grad = rosenbrock(X)
            # no preconditioning in this test: simply copy grad in
            # grad_preco
            grad_preco = np.copy(grad)

    # Helpful console writings
    print('END OF TEST')
    print('FINAL iterate is : ', X)
    print('See the convergence history in iterate_CG.dat')


if __name__ == '__main__':
    main()
