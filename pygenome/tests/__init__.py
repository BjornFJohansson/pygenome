#!/usr/bin/env python
# -*- coding: utf-8 -*-


import sys
import os
import shutil
import logging
import tempfile
import platform
import pytest

# Pytest plugin for measuring coverage.

#- coverage
#- python-coveralls
#- pytest-cov

def run_pytest():
    
    pyversion = platform.python_version()
    
    print('python {}'.format(pyversion))

    args = sys.argv

    if os.getenv("CI"):
        ci = os.getenv("TRAVIS") or os.getenv("APPVEYOR") or os.getenv("CI_NAME") or os.getenv("CIRCLECI")
        print("Tests run on continuous integration server", ci)
        cwd = os.getcwd()
        print("current working directory:", cwd)
        
        os.environ["pygenome_data_dir"]    = os.path.join(cwd,"DATA_"+pyversion)
        os.environ["pygenome_log_dir"]     = os.path.join(cwd,"DATA_"+pyversion)
        os.environ["pygenome_config_dir"]  = os.path.join(cwd,"DATA_"+pyversion)

        # create data directory if not present
        try:
            os.makedirs( os.environ["pygenome_data_dir"] )
        except OSError:
            if os.path.isdir( os.environ["pygenome_data_dir"] ):
                pass
            else:
                raise

        print('os.environ["pygenome_data_dir"] = ',  os.environ["pygenome_data_dir"])
        print('os.environ["pygenome_log_dir"] = ',   os.environ["pygenome_log_dir"])
        print('os.environ["pygenome_config_dir"] = ',os.environ["pygenome_config_dir"])
        print('')

    else:
        print("Tests run locally")
        os.environ["pygenome_data_dir"]    = tempfile.mkdtemp(prefix="pygenome_data_dir_")
        os.environ["pygenome_log_dir"]     = tempfile.mkdtemp(prefix="pygenome_log_dir_")
        os.environ["pygenome_config_dir"]  = tempfile.mkdtemp(prefix="pygenome_config_dir_")
        os.environ["pygenome_loglevel"]    = str( logging.DEBUG )
        import pygenome
        from importlib import reload
        reload(pygenome)
       
    print('tests')
    try:
        import pytest_cov
    except ImportError:
        print("pytest_cov NOT installed")
        args = []
    else:
        del pytest_cov
        args = ["--cov=pygenome", "--cov-report=html", "--cov-report=xml"]
        print("pytest_cov installed")
    try:
        import nbval
    except ImportError:
        print("nbval NOT installed")
    else:
        del nbval
        print("nbval installed")
        args.append("--nbval")

    args = [ "-v", "-s"] + args  # "."
    cwd = os.getcwd()
    my_package_file = os.path.abspath(__file__)
    my_package_dir = os.path.dirname(my_package_file)
    os.chdir(my_package_dir)
    
    #env_name = os.environ["ENV_NAME"]
    #env_no = os.environ["ENV_NUMBER"]
    #pytest_args=(env_name, env_no)
    
    
    
    pytest.main(args)    # cmdline.
    os.chdir(cwd)
    try:
        shutil.copy(os.path.join("tests","coverage.xml"), "coverage.xml")
    except FileNotFoundError:
        pass
    args = ["pygenome", "--doctest-modules", "-v", "-s"]
    print('doctests')
    #pytest.cmdline.main(args)
    print('done!')