#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''This module provides access to the Saccharomyces cerevisiae genome from Python.
   Sequences can be accessed as Bio.SeqRecord objects provided by Biopython.
'''


import os            as _os
import pickle        as _pickle

data_dir = _os.path.join( _os.getenv("pygenome_data_dir"), "Saccharomyces_cerevisiae")

gene = _pickle.load( open( _os.path.join(data_dir, "gene.pickle"), "rb" ) )



if __name__=="__main__":
    pass