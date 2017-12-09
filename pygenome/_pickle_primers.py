#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
gn.JEN1.sgdpage
       .shortdescription
       .longdescription

Standard Name
GDH3
Systematic Name
YAL062W
SGD ID
S000000058
Aliases
FUN51
Feature Type
ORF , Verified
Description
NADP(+)-dependent glutamate dehydrogenase; synthesizes glutamate from ammonia and alpha-ketoglutarate; rate of alpha-ketoglutarate utilization differs from Gdh1p; expression regulated by nitrogen and carbon sources; GDH3 has a paralog, GDH1, that arose from the whole genome duplication 1 2 3 4
Name Description
Glutamate DeHydrogenase
Paralog
GDH1 1
'''

import os            as _os
import pickle        as _pickle
import csv           as _csv
import collections   as _collections

data_dir = _os.path.join( _os.getenv("pygenome_data_dir"), "Saccharomyces_cerevisiae")
from  pygenome._data import _data_files

_primertuple = _collections.namedtuple("_primertuple", '''rec_num
                                                          ORF_name
                                                          deletion_alias
                                                          essential
                                                          A_confirmation_primer_sequence
                                                          B_confirmation_primer_sequence
                                                          C_confirmation_primer_sequence
                                                          D_confirmation_primer_sequence
                                                          UPTAG_primer_sequence
                                                          DNTAG_primer_sequence
                                                          UPstream45_primer_sequence
                                                          DNstream45_primer_sequence
                                                          UPstream90_primer_sequence
                                                          DNstream90_primer_sequence
                                                          AB_wt_PCR
                                                          AkanB_del_PCR
                                                          CD_wt_PCR
                                                          DkanC_del_PCR
                                                          AD_wt
                                                          AD_del
                                                          AB_del_PCR
                                                          CD_del_PCR
                                                          UPTAG_sequence_20mer
                                                          DNTAG_sequence_20mer''')

def pickle_primers():
    _fn = _os.path.join( _os.getenv("pygenome_data_dir"), "Saccharomyces_cerevisiae",_data_files[16])
    with open(_fn, 'rt') as csvfile:
        rd = _csv.reader(csvfile, delimiter='\t')
        field_names = [x.strip() for x in next(rd)]
        next(rd)
        _primers = _collections.defaultdict(tuple)
        for line_ in rd:
            v = _primertuple(*[x.strip() for x in line_])
            _primers[v.ORF_name] = v
    _pickle.dump( _primers, open(_os.path.join(data_dir,"primers.pickle"), "wb" ), -1 )
    
_not_done_tuple=_collections.namedtuple("_not_done_tuple", "ORF_name Gene_name SGD_class")

def pickle_orfs_not_deleted():
    _fn = _fn = _os.path.join( _os.getenv("pygenome_data_dir"), "Saccharomyces_cerevisiae",_data_files[17])
    with open(_fn, 'rt') as csvfile:
        rd = _csv.reader(csvfile, delimiter='\t')
        next(rd)
        next(rd)
        next(rd)
        field_names = [x.strip() for x in next(rd)]
        next(rd)
        _not_done = _collections.defaultdict(tuple)
        for line_ in rd:
            v = _not_done_tuple(*[x.strip() for x in line_])
            _not_done[v.ORF_name] = v
            
    _pickle.dump( _not_done, open(_os.path.join(data_dir,"not_done.pickle"), "wb" ), -1 )
