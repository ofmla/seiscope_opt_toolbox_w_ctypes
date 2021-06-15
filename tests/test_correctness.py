import pathmagic  # noqa
import os
import subprocess
import sys
import filecmp
import numpy as np
import ctypes
from interface import sotb_wrapper, lib_sotb
from ctypes import c_int, c_float, c_bool, POINTER


test_dir = os.path.dirname(os.path.realpath(__file__))
root_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
exe_dir = os.path.join(test_dir, '..', 'bin')
exe_filenames = ['./'+file for file in os.listdir(exe_dir) if file.startswith("test")]


def SafeExtern(mycmd, rundir):
    '''
    Wrapper to call external programs checking the results
    from https://stackoverflow.com/questions/26741316/how-do-i-\
    run-multiple-executables-in-python-one-after-the-other
    '''
    try:  # This allows exceptions to be caught
        retcode = subprocess.call(mycmd, shell=True, cwd=rundir)  # Call ext prog
        if retcode < 0:  # Check the return code errors should be <0
            print(sys.stderr, "Child was terminated by signal", -retcode)
        else:
            print(sys.stderr, "Child returned", retcode)  # For information
    except OSError as e:  # Catch OSErrors and let the user know
        print(sys.stderr, "Execution failed:", e)
        retcode = -1  # Obviously this is an error
    return retcode


def test_run_wrapper_script():
    '''
    Check that the SEISCOPE optimization toolbox (sotb) wrapper work.
    This script demonstrates how to use the wrapper to find the minimum of
    Rosenbrock's banana function (a famous test case for optimization software)
    Basically, the same fortran tests, but only with four methods are executed
    from python. Assertion conditions are used in the script to check that log
    files were created
    '''
    # Create an instance of the SEISCOPE optimization toolbox (sotb) Class.
    sotb = sotb_wrapper()
    floatptr = POINTER(c_float)

    methods = ['PSTD', 'PNLCG', 'LBFGS', 'TRN']
    string = 'See the convergence history in iterate_'

    for i, method in enumerate(methods):
        # Set some fields of the UserDefined derived type in Fortran
        # (ctype structure) - parameter initialization
        n = c_int(2)                     # dimension
        flag = c_int(0)                  # first flag
        sotb.udf.conv = c_float(1e-8)    # tolerance for the stopping criterion
        sotb.udf.print_flag = c_int(1)   # print info in output files
        sotb.udf.debug = c_bool(False)   # level of details for output files
        if method == 'TRN':
            sotb.udf.niter_max = c_int(100)  # maximum iteration number
        else:
            sotb.udf.niter_max = c_int(10000)
        sotb.udf.nls_max = c_int(30)     # max number of linesearch iteration
        sotb.udf.l = c_int(20)
        sotb.udf.niter_max_CG = c_int(5)  # max no. of inner CG iterations

        # intial guess
        X = np.ones(2, dtype=np.float32)*1.5

        # computation of the cost and gradient associated
        # with the initial guess
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
                lib_sotb.rosenbrock(X.ctypes.data_as(floatptr),
                                    ctypes.byref(fcost),
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
            print(string+method[:2]+'.dat')
        elif method == 'PNLCG':
            print(string+method[3:]+'.dat')
        elif method == 'LBFGS':
            print(string+method[1:-1]+'.dat')
        else:
            print(string+method+'.dat and iterate_'+method+'_CG.dat')

    assert os.path.isfile('iterate_ST.dat')
    assert os.path.isfile('iterate_CG.dat')
    assert os.path.isfile('iterate_LB.dat')
    assert os.path.isfile('iterate_TRN.dat')


def test_run_seiscope_opt_tb():
    '''
    Check that the Fortran executables (tests with the Rosenbrock function
    using different gradient-based methods) ran without problems. If at least
    one element of the list is true it means there was some issue
    '''
    assert not any([SafeExtern(exe, exe_dir) for exe in exe_filenames])


def test_byte_by_byte_comp():
    '''
    Test to verify if both set of log files (those obtained by run the Fortran
    tests and by run the script with the wrapper to the library) have the same
    content. A simple byte-by-byte comparison is done with the filecmp library
    http://docs.python.org/library/filecmp.html
    '''
    files1 = [file for file in os.listdir(exe_dir) if file.endswith(".dat")]
    files2 = [file for file in os.listdir(root_dir) if file.endswith(".dat")]
    # just make sure there are some outputs
    assert len(files1) > 0
    assert len(files2) > 0
    files = set(files1).intersection(set(files2))
    cmpfiles = []
    for file in files:
        cmpfiles.append(filecmp.cmp(os.path.join(exe_dir, file),
                                    os.path.join(root_dir, file),
                                    shallow=False))
    assert all(cmpfiles)
