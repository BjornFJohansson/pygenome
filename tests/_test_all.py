#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
test pygenome

https://stackoverflow.com/questions/14270698/get-file-size-using-python-requests-while-only-getting-the-header
'''
import pytest
import os
from pygenome import sg
from pygenome._pretty import pretty_str
from pygenome.systematic_name import _systematic_name
from pygenome._data import  _data_files, _data_urls
from pygenome.update import updater
import requests_mock as rm_module
import shutil, io

data_dir = os.path.join(os.getenv("pygenome_data_dir"), "Saccharomyces_cerevisiae")

@pytest.fixture
def requests_mock(request):
    m = rm_module.Mocker()
    m.start()
    request.addfinalizer(m.stop)
    return m

def test_update(requests_mock):    
    tmp_data_dir = os.path.join(os.getenv("pygenome_data_dir"), "temp")
    shutil.copytree(data_dir, tmp_data_dir)
    
#    for path, url in zip(_data_files,_data_urls): 
#        # set data file to be really old....
#        os.utime(path, times=(path.stat().st_atime, 946684800))  #946684800 = Saturday 1st January 2000 12:00:00 AM
#        #flo = open(path, "rb")
#        flo = io.BytesIO(b"some text data")
#        requests_mock.get(url, 
#                          headers={'last-modified'  : 'Mon, 01 Jan 2001 00:00:00 GMT', #978307200
#                                   'content-length' : "100"}, 
#                          body = flo)
#    updater()  # should "download" the newer files
#    for path, url in zip(_data_files,_data_urls): 
#        # set data file to be the same age remote
#        os.utime(path, times=(path.stat().st_atime, 978307200))
#        flo = io.BytesIO(b"some text datathat will not be used")
#        requests_mock.get(url, 
#                          headers={'last-modified'  : 'Mon, 01 Jan 2001 00:00:00 GMT', #978307200
#                                   'content-length' : "100"}, 
#                          body = flo)
#    updater() # should keep the local files  
    shutil.rmtree(data_dir)
    shutil.copytree(tmp_data_dir, data_dir)