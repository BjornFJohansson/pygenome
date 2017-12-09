#!/usr/bin/env python
# -*- coding: utf-8 -*-

import urllib        as _urllib
import os            as _os
_data_urls  = [ "http://downloads.yeastgenome.org/sequence/S288C_reference/NCBI_genome_source/chr01.gb", 
                "http://downloads.yeastgenome.org/sequence/S288C_reference/NCBI_genome_source/chr02.gb",
                "http://downloads.yeastgenome.org/sequence/S288C_reference/NCBI_genome_source/chr03.gb", 
                "http://downloads.yeastgenome.org/sequence/S288C_reference/NCBI_genome_source/chr04.gb",
                "http://downloads.yeastgenome.org/sequence/S288C_reference/NCBI_genome_source/chr05.gb", 
                "http://downloads.yeastgenome.org/sequence/S288C_reference/NCBI_genome_source/chr06.gb",
                "http://downloads.yeastgenome.org/sequence/S288C_reference/NCBI_genome_source/chr07.gb", 
                "http://downloads.yeastgenome.org/sequence/S288C_reference/NCBI_genome_source/chr08.gb",
                "http://downloads.yeastgenome.org/sequence/S288C_reference/NCBI_genome_source/chr09.gb", 
                "http://downloads.yeastgenome.org/sequence/S288C_reference/NCBI_genome_source/chr10.gb",
                "http://downloads.yeastgenome.org/sequence/S288C_reference/NCBI_genome_source/chr11.gb", 
                "http://downloads.yeastgenome.org/sequence/S288C_reference/NCBI_genome_source/chr12.gb",
                "http://downloads.yeastgenome.org/sequence/S288C_reference/NCBI_genome_source/chr13.gb", 
                "http://downloads.yeastgenome.org/sequence/S288C_reference/NCBI_genome_source/chr14.gb",
                "http://downloads.yeastgenome.org/sequence/S288C_reference/NCBI_genome_source/chr15.gb", 
                "http://downloads.yeastgenome.org/sequence/S288C_reference/NCBI_genome_source/chr16.gb",  
                "http://www-sequence.stanford.edu/group/yeast_deletion_project/Deletion_primers_PCR_sizes.txt", 
                "http://www-sequence.stanford.edu/group/yeast_deletion_project/ORFs_not_available.txt"       ]    

_data_dir = _os.path.join( _os.getenv("pygenome_data_dir"), "Saccharomyces_cerevisiae")

_data_files = [ _urllib.parse.urlparse(url)[2].rpartition("/")[-1] for url in _data_urls]        