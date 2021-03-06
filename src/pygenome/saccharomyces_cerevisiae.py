#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Access to the Saccharomyces cerevisiae genome from Python.

Sequences can be accessed as Bio.SeqRecord objects provided by Biopython.
"""

import sys as _sys
import os as _os

if _sys.version_info >= (3, 8):
    import pickle as _pickle
else:
    import pickle5 as _pickle

_data_dir = _os.path.join(_os.getenv("pygenome_data_dir"),
                          "Saccharomyces_cerevisiae")

sysgenes = _pickle.load(open(_os.path.join(_data_dir, "sysgene.pickle"), "rb"))
stdgenes = _pickle.load(open(_os.path.join(_data_dir, "stdgene.pickle"), "rb"))
