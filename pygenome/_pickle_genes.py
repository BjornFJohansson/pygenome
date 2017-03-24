#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''This module provides access to the Saccharomyces cerevisiae genome from Python.
   Sequences can be accessed as Bio.SeqRecord objects provided by Biopython.
'''
import sys           as _sys
import os            as _os
import pickle        as _pickle

from requests.structures import CaseInsensitiveDict

data_dir = _os.path.join( _os.getenv("pygenome_data_dir"), "Saccharomyces_cerevisiae")

_feature_list = _pickle.load( open(_os.path.join(data_dir, "feature_list.pickle"), "rb" ) )
_standard_to_systematic = _pickle.load( open(_os.path.join(data_dir, "standard_to_systematic.pickle"), "rb" ) )

from pygenome._locus import Gene

def _pickle_genes():   

    gene = CaseInsensitiveDict()
    _sys.stdout.write("\n\nPickle systematic gene names\n")
    for _f in _feature_list:
        gene[_f] = Gene(_f)
        _sys.stdout.write("{} ".format(_f))
    _sys.stdout.write("\n\nPickle standard gene names\n")
    for _f,_g in list(_standard_to_systematic.items()):
        gene[_f] = Gene(_g)
        _sys.stdout.write("{} ".format(_f))
        
    _pickle.dump( gene, open( _os.path.join(data_dir,"gene.pickle"), "wb" ), -1 )

if __name__=="__main__":
    pass