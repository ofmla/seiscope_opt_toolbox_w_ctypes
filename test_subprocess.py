# -*- coding: utf-8 -*-

import subprocess as sp


def readGitVersion():
    try:
        proc = sp.check_output(('git', 'describe', '--long',
                                 '--match', 'v[0-9]*.*'),
                               stderr=sp.STDOUT).splitlines()[0].strip()
    except sp.CalledProcessError as error:
        errorMessage = ">>> Error while executing:\n"\
                       + 'git describe'\
                       + "\n>>> Returned with error:\n"\
                       + str(error.output)
        print("Error: " + errorMessage)

    return proc

if __name__ == '__main__':
    print(readGitVersion())
