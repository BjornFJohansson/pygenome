#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import os
import logging
import tempfile
import platform
import pytest

def main():

    os.environ["pygenome_data_dir"]    = tempfile.mkdtemp(prefix="pygenome_data_dir_")
    os.environ["pygenome_log_dir"]     = tempfile.mkdtemp(prefix="pygenome_log_dir_")
    os.environ["pygenome_config_dir"]  = tempfile.mkdtemp(prefix="pygenome_config_dir_")
    os.environ["pygenome_loglevel"]    = str( logging.DEBUG )

    print("\n\ntests py {}\n\n".format(platform.python_version()))

    args = ["tests",
            "src",
            "--cov=pygenome",
            "--cov-report=html",
            "--cov-report=xml",
            "--import-mode=importlib",
            "--nbval",
            "--current-env",
            "--capture=no",
            "--durations=10",
            "--doctest-modules",
            "-v"]

    result_suite = pytest.cmdline.main(args)

    print("\n\ndone!")

    return result_suite


if __name__ == "__main__":
    result = main()
    sys.exit(result)
