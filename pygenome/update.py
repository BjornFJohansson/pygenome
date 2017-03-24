#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''This module provides access to the Saccharomyces cerevisiae genome from Python.
   Sequences can be accessed as Bio.SeqRecord objects provided by Biopython.
'''

import os            as _os
import urllib        as _urllib
import urllib.parse  as _urlparse
import sys           as _sys
import time          as _time
import datetime      as _datetime

from tqdm import tqdm
import requests

from pydna.readers import read as _read

from pygenome._pFA6a_kanMX4 import plasmid as _plasmid
_pFA6_kanMX4 = _read(_plasmid) # AJ002680

from pygenome._data_files import _chromosome_files, _data_files

data_dir = _os.path.join( _os.getenv("pygenome_data_dir"), "Saccharomyces_cerevisiae")



def _download_chromosome_files(_missing_files=[]):
    '''
    Download the sequence files from Saccharomyces Genome Database (www.sgd.org)
    This is typically only done once.
    '''
    url = _os.getenv("pygenome_base_url")
    _sys.stdout.write("Data files to be downloaded are:\n")

    for _file_ in _missing_files:
        _sys.stdout.write(_file_+"\n")

    _sys.stdout.write("these files will be put in {}\n".format(data_dir))

    for _file_ in sorted(_missing_files):
        _sys.stdout.write("downloading {} ".format(_file_))

        last_modified = _urllib.request.urlopen(_urlparse.urljoin(url, _file_)).headers['last-modified']
        remotedate = _time.strptime(last_modified, '%a, %d %b %Y %H:%M:%S %Z')

        rdate2 = int(_time.mktime(remotedate))
        
        response = requests.get(_urlparse.urljoin(url, _file_), stream=True)
        total = int(response.headers.get('content-length'))
        
        with open(_os.path.join(data_dir ,_file_), 'wb') as f:
            for data in tqdm(response.iter_content(), total=total):
                f.write(data)

        _os.utime(_os.path.join(data_dir ,_file_) ,(rdate2, rdate2))

        _sys.stdout.write("{} successfully downloaded\n".format(_file_))

def _download_data_files(_missing_files=[]):
    '''
    Download data files.
    This is typically only done once.
    '''
    url = _os.getenv("pygenome_primer_url")
    _sys.stdout.write("Data files to be downloaded are:\n")

    for _file_ in _missing_files:
        _sys.stdout.write(_file_+"\n")

    _sys.stdout.write("these files will be put in {}\n".format(data_dir))

    for _file_ in sorted(_missing_files):
    
        from tqdm import tqdm
        import requests
    
        response = requests.get(_urlparse.urljoin(url, _file_), stream=True)
        total = int(response.headers.get('content-length'))
    
        with open(_os.path.join(data_dir ,_file_), 'wb') as f:
            for data in tqdm(response.iter_content(), total=total):
                f.write(data)


def update():
    url = _os.getenv("pygenome_base_url")
    _sys.stdout.write("checking for updated chromosome files at\n{}\n".format(url))
    
    _missing_files = []

    for _file_ in sorted(_chromosome_files.values()):

        last_modified = _urllib.request.urlopen(_urlparse.urljoin(url, _file_)).headers['last-modified']
        remotedate = _time.strptime(last_modified, '%a, %d %b %Y %H:%M:%S %Z')
        # http://stackoverflow.com/questions/5022083/how-can-i-get-the-last-modified-time-with-python3-urllib
        
        if _datetime.datetime(*remotedate[:-2]) > _datetime.datetime.fromtimestamp((_os.path.getmtime(_os.path.join(data_dir, _file_)))):
            _missing_files.append(_file_)
            _sys.stdout.write("{} is available in a newer version\n".format(_file_))
        else:
            _sys.stdout.write("{} is the newest version\n".format(_file_))

    if _missing_files:
        _download_chromosome_files(_missing_files)
    
    url = _os.getenv("pygenome_primer_url")
    _sys.stdout.write("\nchecking for updated primer files at\n{}\n".format(url))
    
    _missing_files = []    
    
    for _file_ in _data_files:
    
        last_modified = _urllib.request.urlopen(_urlparse.urljoin(url, _file_)).headers['last-modified']
        remotedate = _time.strptime(last_modified, '%a, %d %b %Y %H:%M:%S %Z')
    
        if _datetime.datetime(*remotedate[:-2]) > _datetime.datetime.fromtimestamp((_os.path.getmtime(_os.path.join(data_dir, _file_)))):
            _sys.stdout.write("{} is available in a newer version\n".format(_file_))
            _missing_files.append(_file_)
        else:
            _sys.stdout.write("{} is the newest version\n".format(_file_))

    if _missing_files:     
        _download_data_files(_missing_files)

if __name__=="__main__":
    pass