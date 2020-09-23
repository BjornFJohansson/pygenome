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


def main():

    os.environ["pygenome_data_dir"]    = tempfile.mkdtemp(prefix="pygenome_data_dir_")
    os.environ["pygenome_log_dir"]     = tempfile.mkdtemp(prefix="pygenome_log_dir_")
    os.environ["pygenome_config_dir"]  = tempfile.mkdtemp(prefix="pygenome_config_dir_")
    os.environ["pygenome_loglevel"]    = str( logging.DEBUG )

    from pygenome import __file__ as pygenomeinit

    doctestdir = str(pathlib.Path(pygenomeinit).parent)

    args = [ "tests",
             doctestdir,
             "--capture=no",
             "--durations=10",
             "--cov=pygenome",
             "--cov-report=html",
             "--cov-report=xml",
             "--import-mode=importlib",
             "--nbval",
             "--current-env",
             "--doctest-modules",
             "--capture=no",
             "--import-mode=importlib",
             "-vvv"]

    result_suite = pytest.cmdline.main(args)

    return result_suite


if __name__ == "__main__":
    result = main()
    sys.exit(result)
