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

    if "coveralls" in installed:
        print("coveralls-python is installed.")
        args = [
            "--cov=pygenome",
            "--cov-report=html",
            "--cov-report=xml",
            "--import-mode=importlib",
        ]
    else:
        print("coveralls-python NOT installed! (pip install coveralls)")
        args = []

    if "nbval" in installed:
        print("nbval is installed.")
        args.append("--nbval")
        args.append("--current-env")
    else:
        print("nbval NOT installed! (pip install nbval)")

    mainargs = ["tests", "--capture=no", "--durations=10"] + args

    result_suite = pytest.main(mainargs)


    from pygenome import __file__ as pygenomeinit
    doctestdir = str(pathlib.Path(pygenomeinit).parent)
    asciitext("doctests py {}".format(platform.python_version()))

    doctestargs = [
    doctestdir,
    "--doctest-modules",
    "--capture=no",
    "--import-mode=importlib",]

    result_doctest = pytest.main(doctestargs)

    asciitext("done!")

    return result_doctest and result_suite




















    asciitext("tests python {}".format(platform.python_version()))

    installed = {pkg.key for pkg in pkg_resources.working_set}

    if "coveralls" in installed:
        print("coveralls-python is installed.")
        args = [
            "--cov=pydna",
            "--cov-report=html",
            "--cov-report=xml",
            "--import-mode=importlib",
        ]
    else:
        print("coveralls-python NOT installed! (pip install coveralls)")
        args = []

    if "nbval" in installed:
        print("nbval is installed.")
        args.append("--nbval")
        args.append("--current-env")
    else:
        print("nbval NOT installed! (pip install nbval)")


    mainargs = ["tests", "--capture=no", "--durations=10"] + args
    # cwd = os.getcwd()
    # os.chdir("tests")
    result_suite = pytest.cmdline.main(mainargs)
    # os.chdir(cwd)

    from pydna import __file__ as pydnainit

    doctestdir = str(pathlib.Path(pydnainit).parent)

    asciitext("doctests python {}".format(platform.python_version()))
    doctestargs = [
        doctestdir,
        "--doctest-modules",
        "--capture=no",
        "--import-mode=importlib",
    ]
    result_doctest = pytest.cmdline.main(doctestargs)

    asciitext("done!")

    return result_doctest and result_suite


if __name__ == "__main__":
    result = main()
    sys.exit(result)
