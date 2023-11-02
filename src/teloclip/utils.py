#!/usr/bin/env python

import tempfile
import sys
import subprocess
import shutil
import os
from datetime import datetime


def log(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


class Error(Exception):
    pass


def decode(x):
    try:
        s = x.decode()
    except:
        return x
    return s


def _write_script(cmds, script):
    """Write commands into a bash script"""
    f = open(script, "w+")
    for cmd in cmds:
        print(cmd, file=f)
    f.close()


def syscall(cmd, verbose=False):
    """Manage error handling when making syscalls"""
    if verbose:
        print("Running command:", cmd, flush=True)
    try:
        output = subprocess.check_output(cmd, shell=True, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as error:
        print(
            "The following command failed with exit code",
            error.returncode,
            file=sys.stderr,
        )
        print(cmd, file=sys.stderr)
        print("\nThe output was:\n", file=sys.stderr)
        print(error.output.decode(), file=sys.stderr)
        raise Error("Error running command:", cmd)
    if verbose:
        print(decode(output))


def run_cmd(cmds, verbose=False, keeptemp=False):
    """Write and excute script"""
    tmpdir = tempfile.mkdtemp(prefix="tmp.", dir=os.getcwd())
    original_dir = os.getcwd()
    os.chdir(tmpdir)
    script = "run_jobs.sh"
    _write_script(cmds, script)
    syscall("bash " + script, verbose=verbose)
    os.chdir(original_dir)
    if not keeptemp:
        shutil.rmtree(tmpdir)


def getTimestring():
    """Return int only string of current datetime with milliseconds."""
    (dt, micro) = datetime.utcnow().strftime("%Y%m%d%H%M%S.%f").split(".")
    dt = "%s%03d" % (dt, int(micro) / 1000)
    return dt


def dochecks(args):
    """Housekeeping tasks: Create output files/dirs and temp dirs as required."""
    # Make outDir if does not exist else set to current dir.
    if args.temp:
        if not os.path.isdir(args.temp):
            os.makedirs(args.temp)
        tempParent = args.temp
    else:
        tempParent = os.getcwd()
    # Make temp directory
    os.makedirs(os.path.join(tempParent, "temp_" + getTimestring()))
    # Return full path to temp directory
    return tempParent


def missing_tool(tool_name):
    path = shutil.which(tool_name)
    if path is None:
        return [tool_name]
    else:
        return []


def isfile(path):
    """
    Test for existence of input file.
    """
    if not os.path.isfile(path):
        print("Input file not found: %s" % path)
        sys.exit(1)
    else:
        return os.path.abspath(path)
