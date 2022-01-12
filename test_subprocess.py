# -*- coding: utf-8 -*-

import subprocess as sp


def readGitVersion():
    proc = sp.Popen("git describe --long --match v[0-9]*.*",
    				 shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    output, error = proc.communicate()
    if proc.returncode != 0:
        raise Exception("git command call failed %d %s %s" % (process.returncode, output, error))
    return output.splitlines()[0].strip()

if __name__ == '__main__':
    print(readGitVersion())
