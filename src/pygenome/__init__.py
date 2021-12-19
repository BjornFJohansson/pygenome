#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
# pygenome

The pygenome package.

:copyright: Copyright 2013 - 2021 by Björn Johansson. All rights reserved.
:license:   This code is part of the pygenome, distribution and governed by its
            license.  Please see the LICENSE.txt file that should have been
            included as part of this package.

"""
import os as _os
import appdirs as _appdirs
from pygenome._version import version as __version__
import sys as _sys
import subprocess as _subprocess
import prettytable as _prettytable
from pathlib import Path as _Path
from pydna._pretty import pretty_str as _pretty_str
from send2trash import send2trash as _send2trash

__author__ = "Björn Johansson"
__copyright__ = "Copyright 2013-2021, Björn Johansson"
__credits__ = ["Björn Johansson"]
__license__ = "BSD"
__maintainer__ = "Björn Johansson"
__email__ = "bjorn_johansson@bio.uminho.pt"
__status__ = "Development"  # "Production" #"Prototype"
__version__

# create directories

__dirs = [("pygenome_config_dir", _appdirs.user_config_dir),
          ("pygenome_data_dir",   _appdirs.user_data_dir),
          ("pygenome_log_dir",    _appdirs.user_log_dir), ]

for __d, __f in __dirs:
    _os.environ[__d] = _os.getenv(__d, __f("pygenome"))
    _Path(_os.environ[__d]).mkdir(parents=True, exist_ok=True)

# set log level
_os.environ["pygenome_loglevel"] = _os.getenv("pygenome_loglevel",
                                              "DEBUG")


def open_data_dir():
    """docstring."""
    _open_folder(_os.environ["pygenome_data_dir"])


def delete_data_dir():
    """docstring."""
    _send2trash(_os.environ["pygenome_data_dir"])
    _Path(_os.environ["pygenome_data_dir"]).mkdir(parents=True, exist_ok=True)


def open_config_dir():
    """docstring."""
    _open_folder(_os.environ["pygenome_config_dir"])


def open_log_dir():
    """docstring."""
    _open_folder(_os.environ["pygenome_log_dir"])


def _open_folder(pth):
    print(pth)
    if _sys.platform == 'win32':
        _subprocess.Popen(['start', pth], shell=True)
    elif _sys.platform == 'darwin':
        _subprocess.Popen(['open', pth])
    else:
        try:
            _subprocess.Popen(['xdg-open', pth])
        except OSError:
            return "no folder to open."


def get_env():
    """docstring."""
    _table = _prettytable.PrettyTable(["Variable", "Value"])
    _table.set_style(_prettytable.DEFAULT)
    _table.align["Variable"] = "l"  # Left align
    _table.align["Value"] = "l"  # Left align
    _table.padding_width = 1  # One space between column edges and contents
    for d, _ in __dirs:
        _table.add_row((d, _os.environ[d]))
    return _pretty_str(_table)
