# -*- coding: utf-8 -*-

import subprocess as sp


def readGitVersion():
    try:
        proc = sp.check_output(('git', 'describe', '--long',
                                 '--match', 'v[0-9]*.*'),
                               stderr=sp.STDOUT)
    except sp.CalledProcessError as error:
        errorMessage = ">>> Error while executing:\n"\
                       + command\
                       + "\n>>> Returned with error:\n"\
                       + str(error.output)
        print("Error: " + errorMessage)

    return proc.splitlines()[0].strip()

if __name__ == '__main__':
    print(readGitVersion())
