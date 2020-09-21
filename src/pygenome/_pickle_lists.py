#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''This module provides access to the Saccharomyces cerevisiae genome from Python.
   Sequences can be accessed as Bio.SeqRecord objects provided by Biopython.
'''

import platform as _platform

if _platform.python_version().startswith("3.8"):
    import pickle        as _pickle
else:
    import pickle5       as _pickle

import os                               as _os
from Bio             import SeqIO       as _SeqIO
from pygenome._pretty   import pretty_str  as _ps
from  pygenome._data import _data_files

data_dir = _os.path.join( _os.getenv("pygenome_data_dir"), "Saccharomyces_cerevisiae")

def _pickle_lists():
    try:
        _feature_list = _pickle.load( open(_os.path.join(data_dir, "feature_list.pickle"), "rb" ) )
        _standard_to_systematic = _pickle.load( open(_os.path.join(data_dir, "standard_to_systematic.pickle"), "rb" ) )
        _systematic_to_standard = _pickle.load( open(_os.path.join(data_dir, "systematic_to_standard.pickle"), "rb" ) )
        _systematic_to_genbank_accession = _pickle.load( open(_os.path.join(data_dir, "systematic_to_genbank_accession.pickle"), "rb" ) )
        _systematic_to_description = _pickle.load( open(_os.path.join(data_dir, "systematic_to_description.pickle"), "rb" ) )
    except IOError:
        _feature_list  = []
        _standard_to_systematic = {}
        _systematic_to_genbank_accession = {}
        _systematic_to_description = {}

    for _f in _data_files[:16]:
        _krom  =  _SeqIO.read(_os.path.join( _os.getenv("pygenome_data_dir"), "Saccharomyces_cerevisiae", _f), "gb")
        _features = [f for f in _krom.features if f.type=="CDS"]

        _systematic_to_genbank_accession.update( {_ps(f.qualifiers['locus_tag'][0]) : _ps("{} REGION: ".format(_krom.id)+{ 1:"{}..{}".format(f.location.start+1, f.location.end), -1:"complement({}..{})".format(f.location.start+1, f.location.end)}[f.location.strand]) for f in _features })
        _feature_list.extend( [_ps(f.qualifiers['locus_tag'][0]) for f in _features] )
        _standard_to_systematic.update( {_ps(f.qualifiers['gene'][0]):_ps(f.qualifiers['locus_tag'][0]) for f in _features if "gene" in list(f.qualifiers.keys())} )
        _systematic_to_standard = {v: k for k, v in list(_standard_to_systematic.items())}

        for f in _features:
            try:
                description = f.qualifiers["note"][0]
            except KeyError:
                description = f.qualifiers["product"][0]
            _systematic_to_description[_ps(f.qualifiers['locus_tag'][0])] = _ps(description)

        _pickle.dump( _feature_list, open(_os.path.join(data_dir,"feature_list.pickle"), "wb" ), -1 )
        _pickle.dump( _standard_to_systematic, open(_os.path.join(data_dir,"standard_to_systematic.pickle"), "wb" ), -1 )
        _pickle.dump( _systematic_to_genbank_accession, open(_os.path.join(data_dir,"systematic_to_genbank_accession.pickle"), "wb" ), -1 )
        _pickle.dump( _systematic_to_standard, open(_os.path.join(data_dir,"systematic_to_standard.pickle"), "wb" ), -1 )
        _pickle.dump( _systematic_to_description, open(_os.path.join(data_dir,"systematic_to_description.pickle"), "wb" ), -1 )
