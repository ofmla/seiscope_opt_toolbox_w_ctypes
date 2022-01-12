# -*- coding: utf-8 -*-

import subprocess


def readGitVersion():
    proc = subprocess.Popen(('git', 'describe', '--long',
                             '--match', 'v[0-9]*.*'),
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = proc.communicate()
    if proc.returncode != 0:
        raise Exception("git command call failed %d %s %s" % (proc.returncode, output, error))
    return output.splitlines()[0].strip()

if __name__ == '__main__':
    print(readGitVersion())
