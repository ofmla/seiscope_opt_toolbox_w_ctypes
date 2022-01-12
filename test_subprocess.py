# -*- coding: utf-8 -*-

import subprocess as sp


def readGitVersion():
    try:
        proc = sp.check_output(('git', 'describe', '--long',
                               '--match', 'v[0-9]*.*'),
                               stderr=sp.STDOUT).splitlines()[0].strip()
        return proc
    except sp.CalledProcessError as error:
        errorMessage = ">>> Error while executing:\n"\
                       + 'git describe'\
                       + "\n>>> Returned with error:\n"\
                       + str(error.output)
        print("Error: " + errorMessage)
        return None


def git_describe_run():
    proc = sp.run(['git', 'describe', '--long', '--match', 'v[0-9]*.*'],
                  capture_output=True, text=True,
                  check=False)
    print("stdout:", proc.stdout)
    print("stderr:", proc.stderr)

if __name__ == '__main__':
    #print(readGitVersion())
    git_describe_run()
