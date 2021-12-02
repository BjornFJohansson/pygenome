#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""doctsring."""

from Bio.SeqFeature import SeqFeature as Sf
from Bio.SeqFeature import FeatureLocation as Fl

import csv

strand = {True: 1, False: -1}

fss = []
fs = Sf()

chrs = ("chr16.tbl",
        "chr15.tbl",
        "chr14.tbl",
        "chr13.tbl",
        "chr12.tbl",
        "chr11.tbl",
        "chr10.tbl",
        "chr09.tbl",
        "chr08.tbl",
        "chr07.tbl",
        "chr06.tbl",
        "chr05.tbl",
        "chr04.tbl",
        "chr03.tbl",
        "chr02.tbl",
        "chr01.tbl")

for chr in chrs:
    with open(chr) as csvfile:
        csvr = csv.reader(csvfile, delimiter='\t')
        first_row = next(csvr)
        feats = []
        while True:
            try:
                row = next(csvr)
            except StopIteration:
                break
            if row[0]:
                fss.append(fs)
                b, e, *key = row
                b, e = sorted((int(b), int(e)))
                fs = Sf(Fl(b-1, e),
                        type=(key[0:1] or [""])[0],
                        strand=strand[b <= e])
            else:
                _, _, _, qkey, *qval = row
                fs.qualifiers[qkey] = qval
