# -*- coding: utf-8 -*-

import subprocess as sp
import os
import pathlib

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
    proc = sp.run(['git', 'status'],
                  capture_output=True, text=True,
                  check=True, env={"PATH": "/usr/bin"})
    print("stdout:", proc.stdout)
    print("stderr:", proc.stderr)
    print(sp.check_output(["git", "describe", "--always"]).strip().decode())

def get_git_revision(base_path):
    git_dir = pathlib.Path(base_path) / '.git'
    with (git_dir / 'HEAD').open('r') as head:
        ref = head.readline().split(' ')[-1].strip()

    with (git_dir / ref).open('r') as git_hash:
        return git_hash.readline().strip()

if __name__ == '__main__':
    print(readGitVersion())
    git_describe_run()
    print(get_git_revision(os.path.dirname(__file__)))
