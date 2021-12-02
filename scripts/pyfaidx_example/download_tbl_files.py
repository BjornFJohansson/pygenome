#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import urllib as _urllib
import os as _os

data_urls = """\
http://downloads.yeastgenome.org/sequence/S288C_reference/NCBI_genome_source/chr01.tbl
http://downloads.yeastgenome.org/sequence/S288C_reference/NCBI_genome_source/chr02.tbl
http://downloads.yeastgenome.org/sequence/S288C_reference/NCBI_genome_source/chr03.tbl
http://downloads.yeastgenome.org/sequence/S288C_reference/NCBI_genome_source/chr04.tbl
http://downloads.yeastgenome.org/sequence/S288C_reference/NCBI_genome_source/chr05.tbl
http://downloads.yeastgenome.org/sequence/S288C_reference/NCBI_genome_source/chr06.tbl
http://downloads.yeastgenome.org/sequence/S288C_reference/NCBI_genome_source/chr07.tbl
http://downloads.yeastgenome.org/sequence/S288C_reference/NCBI_genome_source/chr08.tbl
http://downloads.yeastgenome.org/sequence/S288C_reference/NCBI_genome_source/chr09.tbl
http://downloads.yeastgenome.org/sequence/S288C_reference/NCBI_genome_source/chr10.tbl
http://downloads.yeastgenome.org/sequence/S288C_reference/NCBI_genome_source/chr11.tbl
http://downloads.yeastgenome.org/sequence/S288C_reference/NCBI_genome_source/chr12.tbl
http://downloads.yeastgenome.org/sequence/S288C_reference/NCBI_genome_source/chr13.tbl
http://downloads.yeastgenome.org/sequence/S288C_reference/NCBI_genome_source/chr14.tbl
http://downloads.yeastgenome.org/sequence/S288C_reference/NCBI_genome_source/chr15.tbl
http://downloads.yeastgenome.org/sequence/S288C_reference/NCBI_genome_source/chr16.tbl
"""
 
data_files = [bn(up(url).path) for url in data_urls.splitlines()] 

import requests

for fn, url in zip(data_files, data_urls.splitlines()):
    r = requests.get(url, allow_redirects=True)
    with open(fn, "wb") as f:
        f.write(r.content)
        
        
        
        