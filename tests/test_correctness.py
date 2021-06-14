import os
import subprocess
import pytest
import runpy
import sys
import filecmp


test_dir = os.path.dirname(os.path.realpath(__file__))
exe_dir = os.path.join(test_dir, '..', 'src', 'bin')
exe_filenames = ['./'+file for file in os.listdir(exe_dir) if file.startswith("test")]
script_dir = os.path.join(test_dir, '..', 'examples')


def SafeExtern(mycmd, rundir):
    '''
    Wrapper to call external programs checking the results
    from https://stackoverflow.com/questions/26741316/how-do-i-\
    run-multiple-executables-in-python-one-after-the-other
    '''
    try:  # This allows exceptions to be caught
        retcode = subprocess.call(mycmd, shell=True, cwd=rundir)  # Ext prog
        if retcode < 0:  # Check the return code errors should be <0
            print(sys.stderr, "Child was terminated by signal", -retcode)
        else:
            print(sys.stderr, "Child returned", retcode)  # For information
    except OSError as e:  # Catch OSErrors and let the user know
        print(sys.stderr, "Execution failed:", e)
        retcode = -1  # Obviously this is an error
    return retcode


def test_run_seiscope_opt_tb():
    '''
    Check that the Fortran executables (tests with the Rosenbrock function
    using different gradient-based methods) ran without problems. If at least
    one element of the list is true it means there was some issue
    '''
    assert not any([SafeExtern(exe, exe_dir) for exe in exe_filenames])


def test_run_wrapper_script():
    '''
    Run a script with assert-statements as pytest test. The script shows how to
    use the wrapper to perform the optimization of the Rosenbrock function.
    Basically, the same fortran tests, but only with four methods are executed
    from python. Assertion conditions are used in the script to check that log
    files were created
    '''
    runpy.run_path(os.path.join(script_dir, 'example.py'))


def test_byte_by_byte_comp():
    '''
    Test to verify if both set of log files (those obtained by run the Fortran
    tests and by run the script with the wrapper to the library) have the same
    content. A simple byte-by-byte comparison is done with the filecmp library
    http://docs.python.org/library/filecmp.html
    '''
    files1 = [file for file in os.listdir(exe_dir) if file.endswith(".dat")]
    files2 = [file for file in os.listdir(script_dir) if file.endswith(".dat")]
    files = set(files1).intersection(set(files2))
    cmpfiles = []
    for file in files:
        cmpfiles.append(filecmp.cmp(os.path.join(exe_dir, file),
                                    os.path.join(script_dir, file),
                                    shallow=False))
    assert all(cmpfiles)

