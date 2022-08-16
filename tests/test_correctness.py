"""Test sotb-wrapper."""
import filecmp
import os
from pathlib import Path
import subprocess
import sys

import numpy as np

from sotb_wrapper import interface

root_dir = Path(__file__).parent.parent
exe_dir = root_dir / "sotb_wrapper"
exe_filenames = ["./" + file for file in os.listdir(exe_dir) if file.startswith("test")]


def SafeExtern(mycmd, rundir):
    """Wrapper to call external programs checking the results.

    from https://stackoverflow.com/questions/26741316/how-do-i-\
    run-multiple-executables-in-python-one-after-the-other

    Args:
        mycmd: External command to be executed
        rundir: Running directory of command

    Returns:
        retcode: Returncode from the executed process
    """
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
    """Check that the SEISCOPE optimization toolbox (sotb) wrapper work.

    This script demonstrates how to use the wrapper to find the minimum of
    Rosenbrock's banana function (a famous test case for optimization software)
    Basically, the same fortran tests, but only with four methods are executed
    from python. Assertion conditions are used in the script to check that log
    files were created
    """
    # Create an instance of the SEISCOPE optimization toolbox (sotb) Class.
    sotb = interface.sotb_wrapper()

    methods = ["PSTD", "PNLCG", "LBFGS", "TRN"]
    string = "See the convergence history in iterate_"
    comments = {
        "PSTD": string + methods[0][1:-1] + ".dat",
        "PNLCG": string + methods[1][3:] + ".dat",
        "LBFGS": string + methods[2][:2] + ".dat",
        "TRN": string + methods[-1] + ".dat and iterate_" + methods[-1] + "_CG.dat",
    }

    def post_comment(comment):
        print(comment)

    for method in methods:
        n = 2  # dimension
        flag = 0  # first flag

        conv = 1e-8  # tolerance for the stopping criterion
        print_flag = 1  # print info in output files
        debug = 0  # level of details for output files
        # maximum iteration number
        if method == "TRN":
            niter_max = 100
        else:
            niter_max = 10000
        nls_max = 20  # max number of linesearch iteration
        l = 20  # maximum number of stored pairs used for the l-BFGS approximation
        niter_max_CG = 5  # maximum number of inner conjugate gradient iterations

        # intial guess
        X = np.ones(n, dtype=np.float32) * 1.5

        # we need to set out another two arrays for TRN method
        d = np.zeros(n, dtype=np.float32)
        Hd = np.zeros(n, dtype=np.float32)

        # computation of the cost and gradient associated
        # with the initial guess
        grad, fcost = sotb.rosenbrock(X)

        # parameter initialization
        # reset udf dict as only a few itens are initialized by 'set_inputs' function
        sotb.udf = interface.UserDefined()
        sotb.set_inputs(
            fcost,
            niter_max,
            nls_max=nls_max,
            conv=conv,
            print_flag=print_flag,
            l=l,
            niter_max_CG=niter_max_CG,
            debug=debug,
        )

        # copy of grad in grad_preco: no preconditioning in this test
        # It is only necessary for PNLCG and PSTD
        # conditional expression is avoided so as to keep down function complexity
        grad_preco = np.zeros(2, dtype=np.float32)
        grad_preco = np.copy(grad)

        # optimization loop: while convergence not reached or
        # linesearch not failed, iterate

        while flag != 2 and flag != 4:
            if method == "PSTD":
                flag = sotb.PSTD(n, X, fcost, grad, grad_preco, flag)
            elif method == "PNLCG":
                flag = sotb.PNLCG(n, X, fcost, grad, grad_preco, flag)
            elif method == "LBFGS":
                flag = sotb.LBFGS(n, X, fcost, grad, flag)
            else:
                flag = sotb.TRN(n, X, fcost, grad, d, Hd, flag)
            if flag == 1:
                # compute cost and gradient at point x
                grad, fcost = sotb.rosenbrock(X)
                # no preconditioning in this test: simply copy grad in grad_preco
                # This is not strictly required for LBFGS and TRN
                # We're only avoiding the use of a conditional statement here
                grad_preco = np.copy(grad)
            elif flag == 7:
                # compute d by the Hessian operator and store in Hd
                Hd = sotb.rosenbrock_hess(X, d)

        # Helpful console writings
        print("END OF TEST")
        print("FINAL iterate is : ", X)
        post_comment(comments[method])

    assert os.path.isfile("iterate_ST.dat")
    assert os.path.isfile("iterate_CG.dat")
    assert os.path.isfile("iterate_LB.dat")
    assert os.path.isfile("iterate_TRN.dat")
    assert os.path.isfile("iterate_TRN_CG.dat")


def test_run_seiscope_opt_tb():
    """Run Fortran execs and check for problems.

    Check that the Fortran executables (tests with the Rosenbrock function
    using different gradient-based methods) ran without problems. If at least
    one element of the list is true it means there was some issue
    """
    assert not any([SafeExtern(exe, exe_dir) for exe in exe_filenames])


def test_byte_by_byte_comp():
    """Compare files.

    Test to verify if both set of log files (those obtained by run the Fortran
    tests and by run the script with the wrapper to the library) have the same
    content. A simple byte-by-byte comparison is done with the filecmp library
    http://docs.python.org/library/filecmp.html
    """
    files1 = [file for file in os.listdir(exe_dir) if file.endswith(".dat")]
    files2 = [file for file in os.listdir(root_dir) if file.endswith(".dat")]
    # just make sure there are some outputs
    assert len(files1) > 0
    assert len(files2) > 0
    files = set(files1).intersection(set(files2))
    cmpfiles = []
    for file in files:
        cmpfiles.append(
            filecmp.cmp(
                os.path.join(exe_dir, file),
                os.path.join(root_dir, file),
                shallow=False,
            )
        )
    assert all(cmpfiles)
