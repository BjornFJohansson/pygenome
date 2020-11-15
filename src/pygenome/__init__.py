#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os               as _os
import sys              as _sys
import subprocess       as _subprocess
import logging          as _logging
import logging.handlers as _handlers
import appdirs          as _appdirs
import configparser     as _configparser
import prettytable      as _prettytable
import shutil           as _shutil
import zipfile          as _zipfile
from pygenome._pretty   import pretty_str  as _pretty_str
from pkg_resources      import resource_filename as _resource_filename

from pygenome import _version



'''
# pygenome

The pygenome package.

:copyright: Copyright 2013 - 2020 by Björn Johansson. All rights reserved.
:license:   This code is part of the pygenome, distribution and governed by its
            license.  Please see the LICENSE.txt file that should have been included
            as part of this package.

'''

__author__       = "Björn Johansson"
__copyright__    = "Copyright 2013-2017, Björn Johansson"
__credits__      = ["Björn Johansson"]
__license__      = "BSD"
__maintainer__   = "Björn Johansson"
__email__        = "bjorn_johansson@bio.uminho.pt"
__status__       = "Development" # "Production" #"Prototype"



__version__ = _version.version
del _version
_sys.modules.pop("pygenome._version", None)

# create config directory
_os.environ["pygenome_config_dir"] = _os.getenv("pygenome_config_dir", _appdirs.user_config_dir("pygenome"))
try:
    _os.makedirs( _os.environ["pygenome_config_dir"] )
except OSError:
    if not _os.path.isdir( _os.environ["pygenome_config_dir"] ):
        raise

# set path for the pygenome.ini file
_ini_path = _os.path.join( _os.environ["pygenome_config_dir"], "pygenome.ini" )

# initiate a config parser instance
_parser = _configparser.ConfigParser()

# if a pygenome.ini exists, it is read
if _os.path.exists(_ini_path):
    _parser.read(_ini_path)
else: # otherwise it is created with default settings
    with open(_ini_path, 'w') as f:
        _parser["main"] = { 'loglevel'  : str(_logging.WARNING),
                            'data_dir'  : _appdirs.user_data_dir("pygenome"),
                            'log_dir'   : _appdirs.user_log_dir("pygenome"),
                            'base_url'  : "http://downloads.yeastgenome.org/sequence/S288C_reference/NCBI_genome_source/",
                            'primer_url': "http://www-sequence.stanford.edu/group/yeast_deletion_project/"}
        _parser.write(f)

# pygenome related environmental variables are set from pygenome.ini if they are not already set
_mainsection = _parser["main"]
_os.environ["pygenome_loglevel"] = _os.getenv("pygenome_loglevel", _mainsection.get("loglevel",str(_logging.WARNING)))
_os.environ["pygenome_data_dir"] = _os.getenv("pygenome_data_dir", _mainsection.get("data_dir",_appdirs.user_data_dir("pygenome")))
_os.environ["pygenome_log_dir"]  = _os.getenv("pygenome_log_dir",  _mainsection.get("log_dir",_appdirs.user_log_dir("pygenome")))
_os.environ["pygenome_base_url"]  = _os.getenv("pygenome_base_url",  _mainsection.get("base_url", "http://downloads.yeastgenome.org/sequence/S288C_reference/NCBI_genome_source/"))
_os.environ["pygenome_primer_url"]  = _os.getenv("pygenome_primer_url",  _mainsection.get("primer_url", "http://www-sequence.stanford.edu/group/yeast_deletion_project/"))

# create log directory if not present
try:
    _os.makedirs( _os.environ["pygenome_log_dir"] )
    _logmsg = "Created log directory {}".format(_os.environ["pygenome_log_dir"])
except OSError:
    if _os.path.isdir( _os.environ["pygenome_log_dir"] ):
        _logmsg = "Log directory {} found.".format(_os.environ["pygenome_log_dir"])
    else:
        raise

# create logger
_logger = _logging.getLogger("pygenome")
_logger.setLevel( int(_os.environ["pygenome_loglevel"]) )
_hdlr = _handlers.RotatingFileHandler(_os.path.join( _os.environ["pygenome_log_dir"] , 'pygenome.log'), mode='a', maxBytes=10*1024*1024, backupCount=10, encoding='utf-8')
_formatter = _logging.Formatter('%(asctime)s %(levelname)s %(funcName)s %(message)s')
_hdlr.setFormatter(_formatter)
_logger.addHandler(_hdlr)
_logger.info(_logmsg)
_logger.info('Assigning environmental variable pygenome_data_dir = %s', _os.environ["pygenome_data_dir"] )

# create data directory if not present
try:
    _os.makedirs( _os.environ["pygenome_data_dir"] )
    _logger.info("Created data directory %s", _os.environ["pygenome_data_dir"] )
except OSError:
    if _os.path.isdir( _os.environ["pygenome_data_dir"] ):
        _logger.info("data directory %s found", _os.environ["pygenome_data_dir"])
    else:
        raise

class _pygenomeWarning(Warning):
    """pygenome warning.

    pygenome uses this warning (or subclasses of it), to make it easy to
    silence all warning messages:

    >>> import warnings
    >>> from pygenome import _pygenomeWarning
    >>> warnings.simplefilter('ignore', _pygenomeWarning)

    Consult the warnings module documentation for more details.
    """
    pass

class _pygenomeDeprecationWarning(_pygenomeWarning):
    """pygenome deprecation warning.

    pygenome uses this warning instead of the built in DeprecationWarning
    since those are ignored by default since Python 2.7.

    To silence all our deprecation warning messages, use:

    >>> import warnings
    >>> from pygenome import _pygenomeDeprecationWarning
    >>> warnings.simplefilter('ignore', _pygenomeDeprecationWarning)

    Code marked as deprecated will be removed in a future version
    of pygenome. This can be discussed in the pygenome google group:
    https://groups.google.com/forum/#!forum/pygenome

    """
    pass

def open_data_folder():
    _open_folder( _os.environ["pygenome_data_dir"] )

def open_config_folder():
    _open_folder( _os.environ["pygenome_config_dir"] )

def open_log_folder():
    _open_folder( _os.environ["pygenome_log_dir"] )

def _open_folder(pth):
    if _sys.platform=='win32':
        _subprocess.Popen(['start', pth], shell=True)
    elif _sys.platform=='darwin':
        _subprocess.Popen(['open', pth])
    else:
        try:
            _subprocess.Popen(['xdg-open', pth])
        except OSError:
            return "no folder to open."

def delete_data_folder():
    import shutil as _shutil
    _shutil.rmtree(_os.environ["pygenome_data_dir"])


def get_env():
    _table = _prettytable.PrettyTable(["Variable", "Value"])
    _table.set_style(_prettytable.DEFAULT)
    _table.align["Variable"] = "l" # Left align
    _table.align["Value"] = "l" # Left align
    _table.padding_width = 1 # One space between column edges and contents
    for k,v in _os.environ.items():
        if k.startswith("pygenome"):
            _table.add_row([k,v])
    return _pretty_str(_table)


data_dir = _os.path.join( _os.getenv("pygenome_data_dir"), "Saccharomyces_cerevisiae")


try:
    _os.makedirs( data_dir )
except OSError:
    if not _os.path.isdir( data_dir ):
       raise

from  pygenome._data import _data_urls
import urllib as _urllib

_missing_files=[]

for _url in _data_urls:
    _file_= _urllib.parse.urlparse(_url)[2].rpartition("/")[-1]
    if not _os.path.exists(_os.path.join(data_dir, _file_)):
        _logger.warning("file %s missing", _file_)
        _missing_files.append(_file_)

import time as _time
import os as _os

if _missing_files:
    _shutil.copy( _resource_filename("pygenome", "Saccharomyces_cerevisiae.zip"), data_dir )
    zf = _zipfile.ZipFile(_os.path.join(data_dir, "Saccharomyces_cerevisiae.zip"), "r")
    for zi in zf.infolist():
        zf.extract(zi, path=data_dir)
        date_time = _time.mktime(zi.date_time + (0, 0, -1))
        _os.utime(_os.path.join(data_dir, zi.filename), (date_time, date_time))
    zf.close()

if not _os.path.exists(_os.path.join(data_dir, "primers.pickle")):
    print("pickle primers start.")
    from pygenome._pickle_primers import pickle_primers
    pickle_primers()
    print("pickle primers done.")
if not _os.path.exists(_os.path.join(data_dir, "not_done.pickle")):
    print("pickle primers.")
    from pygenome._pickle_primers import pickle_orfs_not_deleted
    pickle_orfs_not_deleted()
    print("pickle primers done.")

_do_pickle_lists = False

if not _os.path.exists(_os.path.join(data_dir, "feature_list.pickle")):
    _do_pickle_lists = True
if not _os.path.exists(_os.path.join(data_dir, "standard_to_systematic.pickle")):
    _do_pickle_lists = True
if not _os.path.exists(_os.path.join(data_dir, "systematic_to_standard.pickle")):
    _do_pickle_lists = True
if not _os.path.exists(_os.path.join(data_dir, "systematic_to_genbank_accession.pickle")):
    _do_pickle_lists = True
if not _os.path.exists(_os.path.join(data_dir, "systematic_to_description.pickle")):
    _do_pickle_lists = True



if _do_pickle_lists:
    print("pickle lists start.")
    from pygenome._pickle_lists import _pickle_lists
    _pickle_lists()
    print("pickle lists done.")
if not _os.path.exists(_os.path.join(data_dir,"stdgene.pickle")) or not _os.path.exists(_os.path.join(data_dir,"sysgene.pickle")):
    print("pickle locus list start.")
    from pygenome._pickle_genes import _pickle_genes
    _pickle_genes()
    print("pickle locus list done.")
