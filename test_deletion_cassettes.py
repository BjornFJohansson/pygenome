#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan 29 08:47:17 2022

@author: bjorn
"""

import pickle
from pygenome.saccharomyces_cerevisiae.S288C import data_dir
from pygenome.saccharomyces_cerevisiae.S288C import gene_dicts
from pydna.amplify import pcr
from pydna.primer import Primer
from pydna.assembly import Assembly
from tqdm import tqdm
from Bio import pairwise2
from Bio.pairwise2 import format_alignment

sysgenes, stdgenes = gene_dicts()

# These genes are in the gb files as of 2022-06-08
sysset = set(sysgenes.keys())
assert len(sysgenes) == len(sysset) == 6003

primers = pickle.load(
    (data_dir/"Deletion_primers_PCR_sizes.pickle").open("rb"))

has_primers_set = set(primers.keys())  # 6132
assert len(has_primers_set) == 6132

# ORFs_not_available is a list of genes not deleted at the conclusion of the
# gene deletion project
ndtuple = pickle.load((data_dir/"ORFs_not_available.pickle").open("rb"))
notdeleted_set = set(ndtuple)
assert len(ndtuple) == len(notdeleted_set) == 529

# deleted_in_genome are genes that are:
# 1. present in the latest iteration of the genome annotation
# 2. has an associated primer set
# 3. not in the ORFs_not_available list
deleted_in_genome = sysset & (has_primers_set - notdeleted_set)
assert len(deleted_in_genome) == 5643


# failedset are genes that are:
# 1. present in the latest iteration of the genome annotation
# 2. has an associated primer set
# 3. in the ORFs_not_available list
failedset = sysset & has_primers_set & notdeleted_set

assert failedset == {"YHR119W"}

cas = sysgenes["YHR119W"].deletion_cassettes[0]

loc = sysgenes["YHR119W"].cassette_integration_locus(cas)

# name|43
#      \/
#      /\
#      43|1667bp_PCR_prod|43
#                         \/
#                         /\
#                         43|name

assert sysgenes["YHR119W"].stdname == "SET1"

# Describe deletion:
# Nislow, C., E. Ray, and L. Pillus. 1997. “SET1, a Yeast Member of the
# Trithorax Family, Functions in Transcriptional Silencing and Diverse
# Cellular Processes.” *Molecular Biology of the Cell* 8 (12) (December):
# 2421–2436.

pFA6a_kanMX4 = pickle.load((data_dir/"pFA6a-kanMX4.pickle").open("rb"))

outdated_primers_A = []
outdated_primers_B = []
outdated_primers_C = [] 
outdated_primers_D = []
no_cassette = []
ups45_primer_error = []
no_integration = []

for name in tqdm(deleted_in_genome):

    for p in primers[name]:

        upt = Primer(p[8])
        dnt = Primer(p[9])
        ups45 = Primer(p[10])
        dns45 = Primer(p[11])

        A = Primer(p[4])
        B = Primer(p[5])
        C = Primer(p[6])
        D = Primer(p[7])

        A.id = f"PrimerA_{name}"
        B.id = f"PrimerB_{name}"
        C.id = f"PrimerC_{name}"
        D.id = f"PrimerD_{name}"

        locus = sysgenes[name].locus(upstream=2000,
                                     downstream=2000)

        if not str(A.seq) in locus:
            outdated_primers_A.append(name)

        if not str(B.reverse_complement().seq) in locus:
            outdated_primers_B.append(name)

        if not str(C.seq) in locus:
            outdated_primers_C.append(name)

        if not str(D.reverse_complement().seq) in locus:
            outdated_primers_D.append(name)

        if str(ups45.seq[-18:]) == "GATGTCCACGAGCTCTCT":
            ups45_primer_error.append(name)   
            ups45 = Primer(ups45.seq[:-18]) + "GATGTCCACGAGGTCTCT"
            ups45.id = f"UPTAG_{name}"
            ups45.description = "GATGTCCACGAGCTCTCT -> GATGTCCACGAGGTCTCT"

        if upt and dnt:
            try:
                inner_cassette = pcr(upt, dnt, pFA6a_kanMX4)
            except ValueError:
                inner_cassette = None
            try:
                outer_cassette = pcr(ups45, dns45, inner_cassette)
            except ValueError:
                outer_cassette = None
        else:
            try:
                inner_cassette = pcr(upt, dns45, pFA6a_kanMX4)
            except ValueError:
                inner_cassette = None
            try:
                outer_cassette = pcr(ups45, dns45, inner_cassette)
            except ValueError:
                outer_cassette = None

        if not outer_cassette:
            no_cassette.append(name)

        asm = Assembly((locus, outer_cassette, locus))

        candidates = asm.assemble_linear()

        try:
            kanmx4_gene = candidates[0]
        except IndexError:
            no_integration.append(name)


assert len(outdated_primers_A) == 23
assert len(outdated_primers_B) == 3
assert len(outdated_primers_C) == 139
assert len(outdated_primers_D) == 152
assert len(no_cassette) == 0
assert len(ups45_primer_error) == 318
assert len(no_integration) == 8

print()
for gene in no_integration:

    print(gene)

    cas = sysgenes[gene].deletion_cassettes.pop()

    f = str(cas.forward_primer.seq)
    r = str(cas.reverse_primer.reverse_complement().seq)
    c = str(sysgenes[gene].locus().seq)

    alignment = pairwise2.align.localms(f, c, 2, -1, -.5, -.1)[0]

    print(format_alignment(*alignment))

    alignment = pairwise2.align.localms(r, c, 2, -1, -.5, -.1)[0]

    print(format_alignment(*alignment))

    print("----")

# https://www.biostars.org/p/143799

"""
YCL005W
  1 GCAACTTGTAGGAGGAGAAAGCAG-TATATAACTAGCCGCAATATG
    ||||| |||||||||||||||||| |||||||||||||||||||||
959 GCAAC-TGTAGGAGGAGAAAGCAGGTATATAACTAGCCGCAATATG
  Score=87

   1 TAGGTGATATTGCAATTACTTCTTCTCATGCACTAACAAGTGAAT
     |||||||||||||||||||||||||||||||||||||||||||||
1765 TAGGTGATATTGCAATTACTTCTTCTCATGCACTAACAAGTGAAT
  Score=90

----
YML100W-A
   1 CTT-TCTCTT-GGAACA-A------GAAA---T-A--G-GA--GCAA-T---T-------G---ACA-G------T-TG-T-------CG-----ATG
     ||| |||||| || ||| |      ||||   | |  | ||  |||| |   |       |   ||| |      | || |       ||     |||
1772 CTTGTCTCTTCGG-ACATATTCATGGAAAACTTGACTGCGAATGCAACTACCTCACATACGCCAACAAGCAAGACTATGCTTAAACCCCGGAAAAATG
  Score=74.6

   1 T------AA-ACT----CTT--G-----CT-GTCTG-TT--------TTCAT------CT-G--TGC-AA-GCA-C-A--TC-C-T--GCCA
     |      || |||    |||  |     || ||||  ||        |||||      || |  ||| || ||| | |  || | |  ||||
1744 TACGACAAATACTGCCACTTTAGATGATCTTGTCT-CTTCGGACATATTCATGGAAAACTTGACTGCGAATGCAACTACCTCACATACGCCA
  Score=75.2

----
YCL001W
  1 AAAAAACCTGCCAAGCCCTGCAGAACAATAACAAGCATGTGAATG
    |||||||||||||||||||||||||||||||||||||||||||||
723 AAAAAACCTGCCAAGCCCTGCAGAACAATAACAAGCATGTGAATG
  Score=90

   1 TAACTATGAGAAGGCAGATTCAAGCATATGATAAAATATAGATAT
     ||||||||||||||||||||||| |||||||||||||||||||||
1476 TAACTATGAGAAGGCAGATTCAA-CATATGATAAAATATAGATAT
  Score=87.5

----
YOL141W
  1 GGCCATAAAATTTCGACAACTATAGTGCACATATCTAATACGATG
    |||||||||||||||||||||||||||||||||||||||||||||
959 GGCCATAAAATTTCGACAACTATAGTGCACATATCTAATACGATG
  Score=90

   1 TAATGCTCCTATCGGGATTCGATATGGTTGCCAGCTCGCTATGTG
     |||||||||||||||||||||| ||||||||||||||||||||||
3086 TAATGCTCCTATCGGGATTCGA-ATGGTTGCCAGCTCGCTATGTG
  Score=87.5

----
YGR216C
  1 CAGGAACCTTCGAATTCAGGG-CGGTTT-ATACAGTTAAATGACATG
    ||||||||||||||||||||| |||||| ||||||||||||||||||
957 CAGGAACCTTCGAATTCAGGGGCGGTTTTATACAGTTAAATGACATG
  Score=89

   1 TAACGTTTTCTTTTATAAGTAAAGTAGCGTTTTCTTGCTTCAGCT
     |||||||||||||||||||||||||||||||||||||||||||||
2828 TAACGTTTTCTTTTATAAGTAAAGTAGCGTTTTCTTGCTTCAGCT
  Score=90

----
YCR026C
  1 TACATATATCTGAAAAGAGGTACAGTATAACACTCTCTATAGATG
    |||||||||||||||||||||||||||||||||||||||||||||
959 TACATATATCTGAAAAGAGGTACAGTATAACACTCTCTATAGATG
  Score=90

   1 TAAGCTTAAGGTGGGACAACT-GCGTTATGAAATAAAAGACGCATC
     ||||||||||||||||||||  ||||||||||||||||||||||||
3227 TAAGCTTAAGGTGGGACAAC-CGCGTTATGAAATAAAAGACGCATC
  Score=87

----
YBL066C
  1 CGTATACAGATTATATTGGCTCTGCGTATACGCATTCTCGTCATG
    |||||||||||||||||||||||||||||||||||||||||||||
959 CGTATACAGATTATATTGGCTCTGCGTATACGCATTCTCGTCATG
  Score=90

   1 TAAAAG-ACGGG-ATATCCACCTCTGAG-TTG--TTCCAAA-GT--TGATATC
     |||||| || || ||||||||||||||| |||  ||||||| ||  |||||||
4183 TAAAAGCAC-GGCATATCCACCTCTGAGCTTGTTTTCCAAAAGTAATGATATC
  Score=84.3

----
YLR451W
  1 AAAAATCGCTTCGTAACATTAATACAAATTCTTTTTGCAATTATG
    |||||||||||||||||||||||||||||||||||||||||||||
959 AAAAATCGCTTCGTAACATTAATACAAATTCTTTTTGCAATTATG
  Score=90

   1 TAAAGTCCTTTTCTTTTTTTGCCGTAATGTTTACTTACCCTCGAA
     ||||||||||||||||||||| |||||||||||||||||||||||
3659 TAAAGTCCTTTTCTTTTTTTG-CGTAATGTTTACTTACCCTCGAA
  Score=87.5

----
"""


