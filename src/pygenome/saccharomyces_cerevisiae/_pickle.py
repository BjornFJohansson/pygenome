#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Access to the Saccharomyces cerevisiae genome from Python.

Sequences can be accessed as Bio.SeqRecord objects provided by Biopython.
"""

import sys as _sys
if _sys.version_info >= (3, 8):
    import pickle
else:
    import pickle5 as pickle
import os
import csv
import collections
from pathlib import Path
from collections import defaultdict as dd

from Bio import SeqIO
from pydna._pretty import pretty_str as ps
from requests.structures import CaseInsensitiveDict
# from pygenome.locus import Gene

data_dir = Path(os.getenv("pygenome_data_dir"))/"Saccharomyces_cerevisiae"

with (data_dir/"feature_list.pickle").open('rb') as fp:
    feature_list = pickle.load(fp)


def _pickle_lists():

    # set the ch_urls and ch_file_names variables
    ch_urls = ""
    ch_file_names = []
    with open(data_dir/"settings_Saccharomyces_cerevisiae.py", "r") as f:
        exec(f.read(), globals(), locals())
    assert ch_urls
    assert ch_file_names

    {'misc_RNA', 'centromere', 'telomere', 'source',
     'mobile_element', 'tRNA', 'misc_feature',
     'rep_origin', 'rRNA', 'gene', 'repeat_region',
     'CDS', 'ncRNA', 'mRNA', 'regulatory'}

    feature_list = []
    standard_to_systematic = {}
    systematic_to_genbank_accession = {}
    systematic_to_description = {}
    features = dd(list)
    
    ac_dct = { 1:"{}..{}", 
              -1:"complement({}..{})"}
    
    for ch in ch_file_names:
        krom = SeqIO.read(data_dir/ch, "gb")
        for f in krom.features:
            features[f.type].append(f)
        
        for CDS in features["CDS"]:
            
            systematic_name = ps(f.qualifiers['locus_tag'][0])
            strand = CDS.location.strand
            genbank_accession = ac_dct[strand].format(f.location.start+1, 
                                                      f.location.end)                  
            standard_name
            


        standard_to_systematic.update({ps(f.qualifiers['gene'][0]):ps(f.qualifiers['locus_tag'][0]) for f in CDS if "gene" in list(f.qualifiers.keys())} )

        systematic_to_standard = {v: k for k, v in list(standard_to_systematic.items())}

        for f in CDS:
            try:
                description = f.qualifiers["note"][0]
            except KeyError:
                description = f.qualifiers["product"][0]
            systematic_to_description[ps(f.qualifiers['locus_tag'][0])] = ps(description)

        pickle.dump(feature_list, open(os.path.join(data_dir,"feature_list.pickle"), "wb" ), pickle.HIGHEST_PROTOCOL )
        pickle.dump(standard_to_systematic, open(os.path.join(data_dir,"standard_to_systematic.pickle"), "wb" ), pickle.HIGHEST_PROTOCOL )
        pickle.dump(systematic_to_genbank_accession, open(os.path.join(data_dir,"systematic_to_genbank_accession.pickle"), "wb" ), pickle.HIGHEST_PROTOCOL )
        pickle.dump(systematic_to_standard, open(os.path.join(data_dir,"systematic_to_standard.pickle"), "wb" ), pickle.HIGHEST_PROTOCOL )
        pickle.dump(systematic_to_description, open(os.path.join(data_dir,"systematic_to_description.pickle"), "wb" ), pickle.HIGHEST_PROTOCOL )


def _pickle_genes():
    
    stdgene = CaseInsensitiveDict()
    sysgene = CaseInsensitiveDict()
    
    data_dir = Path(os.getenv("pygenome_data_dir"))/"Saccharomyces_cerevisiae"
    
    feature_list = pickle.load(Path(data_dir)/"feature_list.pickle".open("rb"))
    standard_to_systematic = pickle.load(Path(data_dir)/"standard_to_systematic.pickle".open("rb"))    

    for f in feature_list:
        sysgene[f] = Gene(f)

    for f, g in list(standard_to_systematic.items()):
        stdgene[f] = Gene(g)

    pickle.dump(sysgene, Path(data_dir)/"sysgene.pickle".open("wb"), pickle.HIGHEST_PROTOCOL)
    pickle.dump(stdgene, Path(data_dir)/"stdgene.pickle".open("wb"), pickle.HIGHEST_PROTOCOL)
    
     
    
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





def pickle_primers():
    """docstring."""
    primertuple = collections.namedtuple("_primertuple",
                                           '''rec_num
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
    

    
    fn = data_dir/"Deletion_primers_PCR_sizes.txt"
    primers = collections.defaultdict(list)
    with open(fn, 'rt') as csvfile:
        rd = csv.reader(csvfile, delimiter='\t')
        # field_names = [x.strip() for x in next(rd)]
        next(rd)
        # results = []
        for line_ in rd:
            v = primertuple(*[x.strip() for x in line_])
            primers[v.ORF_name].append(v)
    pickle.dump(primers, open(data_dir/"primers.pickle", "wb"), pickle.HIGHEST_PROTOCOL)


def pickle_orfs_not_deleted():
    """docstring."""

    not_done_tuple = collections.namedtuple("_not_done_tuple",
                                              "ORF_name Gene_name SGD_class")

    fn = data_dir/"ORFs_not_available.txt"
    with open(fn, 'rt') as csvfile:
        rd = csv.reader(csvfile, delimiter='\t')
        next(rd)
        next(rd)
        next(rd)
        # field_names = [x.strip() for x in next(rd)]
        next(rd)
        not_done = collections.defaultdict(tuple)
        for line_ in rd:
            v = not_done_tuple(*[x.strip() for x in line_])
            not_done[v.ORF_name] = v
    pickle.dump(not_done, open(data_dir/"not_done.pickle", "wb"), pickle.HIGHEST_PROTOCOL)