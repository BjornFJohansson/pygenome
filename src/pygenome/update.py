#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''This module provides access to the Saccharomyces cerevisiae genome from Python.
   Sequences can be accessed as Bio.SeqRecord objects provided by Biopython.
'''

import os            as _os
import urllib        as _urllib
#import sys           as _sys
import calendar      as _calendar
import datetime      as _datetime
from   email.utils   import parsedate_to_datetime
import logging       as _logging
_module_logger = _logging.getLogger("pygenome."+__name__)
from tqdm import tqdm
import requests
import pathlib

#from pygenome._pFA6a_kanMX4 import plasmid as _plasmid
#_pFA6_kanMX4 = _read(_plasmid) # AJ002680

from  pygenome._data import _data_urls, _data_files

data_dir = _os.path.join( _os.getenv("pygenome_data_dir"), "Saccharomyces_cerevisiae")

def updater():
    print("checking online for updated data files.")

    for url,fn in zip( _data_urls, _data_files):

        path = pathlib.Path(data_dir).joinpath(fn)

        try:
            local_last_mod = _datetime.datetime.fromtimestamp(path.stat().st_mtime, tz=_datetime.timezone.utc)
        except FileNotFoundError:
            local_last_mod = _datetime.datetime.fromtimestamp(0, tz=_datetime.timezone.utc)

        response = requests.get(url, stream=True)

        remote_last_mod = parsedate_to_datetime(response.headers.get('last-modified'))

        if local_last_mod > remote_last_mod: # local file is newer! Should probably not happen!
            _module_logger.critical("local file %s %s is newer than remote file %s %s", fn, local_last_mod, url, remote_last_mod )

        if remote_last_mod > local_last_mod:
            print("{} is available in a newer version --downloading".format(fn))

            remote_last_mod_time_stamp = _calendar.timegm(remote_last_mod.timetuple())

            total = int(response.headers.get('content-length'))

            with open(str(path), 'wb') as f:
                for data in tqdm(response.iter_content(), total=total):
                    f.write(data)

            _os.utime(str(path), times=(remote_last_mod_time_stamp,)*2)

            print("{} successfully downloaded".format(fn))
        else:
            print("{} is the newest version ({})".format(fn, local_last_mod.isoformat()))

if __name__=="__main__": # pragma: no cover
    updater()
