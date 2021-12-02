#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from urllib.parse import urlparse as up
from os.path import basename as bn
import os as _os

data_urls = """\
http://downloads.yeastgenome.org/sequence/S288C_reference/NCBI_genome_source/chr01.fsa
http://downloads.yeastgenome.org/sequence/S288C_reference/NCBI_genome_source/chr02.fsa
http://downloads.yeastgenome.org/sequence/S288C_reference/NCBI_genome_source/chr03.fsa
http://downloads.yeastgenome.org/sequence/S288C_reference/NCBI_genome_source/chr04.fsa
http://downloads.yeastgenome.org/sequence/S288C_reference/NCBI_genome_source/chr05.fsa
http://downloads.yeastgenome.org/sequence/S288C_reference/NCBI_genome_source/chr06.fsa
http://downloads.yeastgenome.org/sequence/S288C_reference/NCBI_genome_source/chr07.fsa
http://downloads.yeastgenome.org/sequence/S288C_reference/NCBI_genome_source/chr08.fsa
http://downloads.yeastgenome.org/sequence/S288C_reference/NCBI_genome_source/chr09.fsa
http://downloads.yeastgenome.org/sequence/S288C_reference/NCBI_genome_source/chr10.fsa
http://downloads.yeastgenome.org/sequence/S288C_reference/NCBI_genome_source/chr11.fsa
http://downloads.yeastgenome.org/sequence/S288C_reference/NCBI_genome_source/chr12.fsa
http://downloads.yeastgenome.org/sequence/S288C_reference/NCBI_genome_source/chr13.fsa
http://downloads.yeastgenome.org/sequence/S288C_reference/NCBI_genome_source/chr14.fsa
http://downloads.yeastgenome.org/sequence/S288C_reference/NCBI_genome_source/chr15.fsa
http://downloads.yeastgenome.org/sequence/S288C_reference/NCBI_genome_source/chr16.fsa
"""

data_files = [bn(up(url).path) for url in data_urls.splitlines()] 

import requests

for fn, url in zip(data_files, data_urls.splitlines()):
    r = requests.get(url, allow_redirects=True)
    with open(fn, "wb") as f:
        f.write(r.content)
        
        
        
        