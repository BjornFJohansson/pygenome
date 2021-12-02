#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os as _os
from Bio import SeqIO as _SeqIO
from pathlib import Path as _Path
from pydna.utils import flatten

# to set env variables
import pygenome
pygenome

# directory to find data files
data_dir = _Path(_os.getenv("pygenome_data_dir"))/"Saccharomyces_cerevisiae"

ch_urls = ""
ch_file_names = []

with open(data_dir/"settings_Saccharomyces_cerevisiae.py", "r") as f:
    exec(f.read(), globals(), locals())

assert ch_urls
assert ch_file_names

# all features on all chromosomes combined
features = []
for ch_file_name in ch_file_names:
    krom = _SeqIO.read(data_dir/ch_file_name, "gb")
    features.extend(krom.features)


feature_types = set([f.type for f in features])
print(feature_types)

for feature_type in feature_types:
    print(feature_type, len([f for f in features if f.type == feature_type]))

"""
source 16
                    gene 6418
ncRNA 106
misc_feature 14
                    CDS 6003
tRNA 275
                    mRNA 6009
telomere 32
mobile_element 50
rep_origin 352
repeat_region 384
regulatory 12
misc_RNA 10
rRNA 12
centromere 64
"""

print("gene")
genef = [f for f in features if f.type == "gene"]
assert len(genef) == 6418
qualkeys = set(flatten(f.qualifiers.keys() for f in genef))
for key in sorted(qualkeys):
    y = len([f for f in genef if key in f.qualifiers])
    n = len([f for f in genef if key not in f.qualifiers])
    print(key, y, n)
print()
    
genefdict = {f.qualifiers["locus_tag"][0]: f.qualifiers for f in genef}

[f for f in genefdict if "pseudo" in genefdict[f]]


print("mRNA")
mRNAf = [f for f in features if f.type == "mRNA"]
assert len(mRNAf) == 6009
qualkeys = set(flatten(f.qualifiers.keys() for f in mRNAf))
for key in sorted(qualkeys):
    y = len([f for f in mRNAf if key in f.qualifiers])
    n = len([f for f in mRNAf if key not in f.qualifiers])
    print(key, y, n)
print()

mRNAfdict = {f.qualifiers["locus_tag"][0]: f.qualifiers for f in mRNAf}

[f for f in mRNAfdict if "pseudo" in mRNAfdict[f]]

print("CDS")
cdsf = [f for f in features if f.type == "CDS"]
assert len(cdsf) == 6003
qualkeys = set(flatten(f.qualifiers.keys() for f in cdsf))
for key in sorted(qualkeys):
    y = len([f for f in cdsf if key in f.qualifiers])
    n = len([f for f in cdsf if key not in f.qualifiers])
    print(key, y, n)









"""
gene
db_xref         12 6406
experiment 6 6412
gene 5385 1033
gene_synonym 2017 4401
                            locus_tag 6418 0
note 12 6406
pseudo 18 6400

mRNA
db_xref         12 5997
gene 5209 800
gene_synonym 2014 3995
                           locus_tag 6009 0
product 6009 0
pseudo 12 5997

CDS
EC_number 1415 4588
                            codon_start 6003 0
                            db_xref 6003 0
experiment 5434 569
gene 5209 794
gene_synonym 2012 3991
                            locus_tag 6003 0
note 5940 63
                            product 6003 0
protein_id 5984 19
pseudo 6 5997
ribosomal_slippage 47 5956
translation 5997 6
"""



"""
x[3].__dict__{'location': FeatureLocation(ExactPosition(1806), ExactPosition(2169), strand=-1),
 'type': 'gene',
 'id': '<unknown id>',
 'qualifiers': OrderedDict([('gene', ['PAU8']),
                            ('locus_tag', ['YAL068C'])])}

x[4].__dict__{'location': FeatureLocation(ExactPosition(1806), ExactPosition(2169), strand=-1),
'type': 'mRNA',
'id': '<unknown id>',
'qualifiers': OrderedDict([('gene', ['PAU8']),
                           ('locus_tag', ['YAL068C']),
                           ('product', ['seripauperin PAU8'])])}


x[5].__dict__{'location': FeatureLocation(ExactPosition(1806), ExactPosition(2169), strand=-1),

'type': 'CDS',
'id': '<unknown id>',
'qualifiers': OrderedDict([('gene', ['PAU8']),
('locus_tag', ['YAL068C']),
('experiment',['EXISTENCE:mutant phenotype:GO:0030437 ascospore formation [PMID:12586695]',
  'EXISTENCE:mutant phenotype:GO:0045944 positive regulation of transcription by RNA polymerase II [PMID:12586695]']),              
('note', ['hypothetical protein; member of the seripauperin multigene family encoded mainly in subtelomeric regions']),              
('codon_start', ['1']),              
('product', ['seripauperin PAU8']),              
('protein_id', ['DAA06918.1']),              
('db_xref', ['SGD:S000002142']),              
('translation',               ['MVKLTSIAAGVAAIAATASATTTLAQSDERVNLVELGVYVSDIRAHLAQYYMFQAAHPTETYPVEVAEAVFNYGDFTTMLTGIAPDQVTRMITGVPWYSSRLKPAISSALSKDGIYTIAN'])])}
"""

assert [f.qualifiers["locus_tag"] for f in cdsf+genef]

    
    
