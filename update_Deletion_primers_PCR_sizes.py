#!/usr/bin/env python
# coding: utf-8

from tqdm import tqdm
import pandas as pd

df = pd.read_csv("data/Deletion_primers_PCR_sizes.txt", delimiter="\t", header=0)

df = df.rename(columns=lambda x: x.strip())

df = df.iloc[1:].reset_index(drop=True)

df = df.map(lambda x: x.strip() if isinstance(x, str) else x)

U1 = "GATGTCCACGAGGTCTCT"
U2 = "CGTACGCTGCAGGTCGAC"
D1 = "CGGTGTCGGTCTCGTAG"
D2 = "ATCGATGAATTCGAGCTCG"

mutated_U1 = "GATGTCCACGAGCTCTCT"

number_of_rows = len(df)
assert number_of_rows == 6363

no_dntag = 0
for index, row in df.iloc[1:].iterrows():
    row_dict = row.to_dict()
    uptag = row_dict["UPTAG_primer_sequence"]
    dntag = row_dict["DNTAG_primer_sequence"]
    tag1 = row_dict["UPTAG_sequence_20mer"]
    assert U1 + tag1 + U2 in uptag
    if dntag:
        tag2 = row_dict["DNTAG_sequence_20mer"]
        assert D1 + tag2 + D2 in dntag
    else:
        no_dntag += 1


# 192 of the 6363 rows (~3%) contain no DNTAG_primer_sequence which is also confirmed on the original [website](http://chemogenomics.pharmacy.ubc.ca/GGCN_Lab/SGDP/group/yeast_deletion_project/project_desc.html).

assert no_dntag == 192

from pydna.primer import Primer
from pydna.amplify import pcr
from pydna.readers import read

pFA6 = read("data/pFA6a-kanMX4.gb")

assert pFA6.seguid() == 'cdseguid=FkTDjvMc5QaNMAE_dhD52VXLa6c'

uptag_and_dntag = []
uptag_and_dns45 = []

outer_cassettes = {}

for index, row in tqdm(df.iterrows()):
    row_dict = row.to_dict()

    genename = row_dict["ORF_name"]

    if genename in ["YHR119W",]:
        continue

    upt = Primer(row_dict["UPTAG_primer_sequence"])
    dnt = Primer(row_dict["DNTAG_primer_sequence"])
    ups45 = Primer(row_dict["UPstream45_primer_sequence"])
    dns45 = Primer(row_dict["DNstream45_primer_sequence"])
    ups90 = Primer(row_dict["UPstream90_primer_sequence"])
    dns90 = Primer(row_dict["DNstream90_primer_sequence"])

    upt.id = f"UPTAG_{genename}"
    dnt.id = f"DNTAG_{genename}"
    ups45.id = f"UPstream45_{genename}"
    dns45.id = f"DNstream45_{genename}"
    ups90.id = f"UPstream90_{genename}"
    dns90.id = f"DNstream90_{genename}"

    if upt and dnt:
        inner_cassette = pcr(upt, dnt, pFA6)
        uptag_and_dntag.append(genename)
        strategy = "Strategy 1: UPTAG, DNTAG"
        inner_cassette.annotations["comment"] = strategy
        try:
            outer_cassette = pcr(ups45, dns45, inner_cassette)
        except ValueError:
            assert mutated_U1 in ups45
            inner_cassette.seq = inner_cassette.seq.replace(U1, mutated_U1)
            outer_cassette = pcr(ups45, dns45, inner_cassette)
    elif upt and dns45:
        outer_cassette = pcr(upt, dns45, pFA6)
        uptag_and_dns45.append(genename)
        strategy = "Strategy 2: UPTAG, DNstream"
        outer_cassette.annotations["comment"] = strategy
    if ups90 and dns90:
        try:
            long_cassette = pcr(ups90, dns90, outer_cassette)
        except ValueError:
            print(genename)
    outer_cassettes[genename] = outer_cassette
    break

    x = outer_cassettes["YGR194C"]
    x.locus = "ygr194c::KanMX4"

    print(x.format())
