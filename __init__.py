#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
# pygenome

The pygenome package.

:copyright: Copyright 2013 - 2022 by Björn Johansson. All rights reserved.
:license:   This code is part of the pygenome, distribution and governed by its
            license.  Please see the LICENSE.txt file that should have been
            included as part of this package.

"""
import os as _os
import appdirs as _appdirs
from pygenome._version import version as __version__
import sys as _sys
import subprocess as _subprocess
from pathlib import Path as _Path
from send2trash import send2trash as _send2trash

__author__ = "Björn Johansson"
__copyright__ = "Copyright 2013-2022, Björn Johansson"
__credits__ = ["Björn Johansson"]
__license__ = "BSD"
__maintainer__ = "Björn Johansson"
__email__ = "bjorn_johansson@bio.uminho.pt"
__status__ = "Development"  # "Production" #"Prototype"
__version__

# create directories

pygenome_data_dir = _Path(_appdirs.user_data_dir("pygenome"))

pygenome_data_dir.mkdir(parents=True, exist_ok=True)

# set log level
_os.environ["pygenome_loglevel"] = _os.getenv("pygenome_loglevel",
                                              "DEBUG")


def delete_data_dir():
    """docstring."""
    _send2trash(pygenome_data_dir)
    _Path(pygenome_data_dir).mkdir(parents=True, exist_ok=True)


def open_data_dir():
    """docstring."""
    print(pygenome_data_dir)
    if _sys.platform == 'win32':
        _subprocess.Popen(['start', pygenome_data_dir], shell=True)
    elif _sys.platform == 'darwin':
        _subprocess.Popen(['open', pygenome_data_dir])
    else:
        try:
            _subprocess.Popen(['xdg-open', pygenome_data_dir])
        except OSError:
            return "no folder to open."
