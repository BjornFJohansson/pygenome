#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''This module provides access to the Saccharomyces cerevisiae genome from Python.
   Sequences can be accessed as Bio.SeqRecord objects provided by Biopython.
'''
import sys as _sys
if _sys.version_info >= (3,8):
    import pickle  as _pickle
else:
    import pickle5 as _pickle

import sys           as _sys
import os            as _os

from requests.structures import CaseInsensitiveDict
from pygenome._pickle_primers import _primertuple

stdgene  = CaseInsensitiveDict()
sysgene  = CaseInsensitiveDict()

data_dir = _os.path.join( _os.getenv("pygenome_data_dir"), "Saccharomyces_cerevisiae")

_feature_list = _pickle.load( open(_os.path.join(data_dir, "feature_list.pickle"), "rb" ) )
_standard_to_systematic = _pickle.load( open(_os.path.join(data_dir, "standard_to_systematic.pickle"), "rb" ) )

from pygenome.locus import Gene

def _pickle_genes():

    for _f in _feature_list:
        sysgene[_f] = Gene(_f)

    print("Pickle standard gene names.")
    for _f,_g in list(_standard_to_systematic.items()):
        stdgene[_f] = Gene(_g)

    _pickle.dump( sysgene, open( _os.path.join(data_dir,"sysgene.pickle"), "wb" ), -1 )
    _pickle.dump( stdgene, open( _os.path.join(data_dir,"stdgene.pickle"), "wb" ), -1 )
