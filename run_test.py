#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import os
import logging
import tempfile
import platform
import pytest
import pathlib
import pkg_resources


try:
    from pyfiglet import Figlet
except ImportError:
    asciitext = print
else:

    def asciitext(*args, **kwargs):
        f = Figlet(font="doom")
        print(f.renderText(" ".join(args)), **kwargs)


def main():

    os.environ["pygenome_data_dir"]    = tempfile.mkdtemp(prefix="pygenome_data_dir_")
    os.environ["pygenome_log_dir"]     = tempfile.mkdtemp(prefix="pygenome_log_dir_")
    os.environ["pygenome_config_dir"]  = tempfile.mkdtemp(prefix="pygenome_config_dir_")
    os.environ["pygenome_loglevel"]    = str( logging.DEBUG )

    asciitext("tests py {}".format(platform.python_version()))

    installed = {pkg.key for pkg in pkg_resources.working_set}

    if "pytest-cov" in installed:
        print("pytest-cov is installed.")
        args = [
            "--cov=pygenome",
            "--cov-report=html",
            "--cov-report=xml",
            "--import-mode=importlib",
        ]
    else:
        print("pytest-cov NOT installed! (pip install pytest-cov)")
        args = []

    if "nbval" in installed:
        print("nbval is installed.")
        args.append("--nbval")
        args.append("--current-env")
    else:
        print("nbval NOT installed! (pip install nbval)")

    mainargs = [
                 "--capture=no",
                 "--durations=10",
                 "--doctest-modules",
                 "-v"]

    result_suite = pytest.cmdline.main(mainargs + args)

    asciitext("done!")

    return result_suite


if __name__ == "__main__":
    result = main()
    sys.exit(result)
