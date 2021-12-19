#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Access to the Saccharomyces cerevisiae genome from Python.

Sequences can be accessed as Bio.SeqRecord objects provided by Biopython.
"""
import os
import calendar
from datetime import datetime as dt
from datetime import timezone as tz
from email.utils import parsedate_to_datetime
from tqdm import tqdm
import requests
from pathlib import Path
import logging
import sys

import csv
import collections
from collections import defaultdict as dd
from urllib.parse import urlparse as up
from os.path import basename as bn

from Bio.Seq import reverse_complement as rc

from Bio import SeqIO
from pydna._pretty import pretty_str as ps
from requests.structures import CaseInsensitiveDict
from pyfaidx import Fasta

if sys.version_info >= (3, 8):
    import pickle
else:
    import pickle5 as pickle

# set the ch_urls and ch_file_names variables
chromosome_urls = ""
with open(Path(os.getenv("pygenome_config_dir"))/"chromosome_urls.py", "r") as f:
    exec(f.read(), globals(), locals())
assert chromosome_urls

data_dir = Path(os.getenv("pygenome_data_dir"))/"S288C"

ac_dct = {1: "{}..{}",
          -1: "complement({}..{})"}

# url = chromosome_urls.splitlines()[0]


class Gene():
    """docstring."""

    def __init__(self,
                 sysname,
                 stdname,
                 start,
                 end,
                 strand,
                 accession,
                 pth):
        self.sysname = sysname
        self.stdname = stdname
        self.start = start
        self.end = end
        self.strand = strand
        self.accession = accession
        self.pth = pth
        self.pred = None
        self.succ = None

    def cds(self):
        """docstring."""
        s = Fasta(str(self.pth))[self.accession][self.start:self.end]
        return s if self.strand == 1 else -s

    def upstream_igr(self):
        """docstring."""
        if self.strand == 1:
            s = Fasta(str(self.pth))[self.accession][self.pred.end,
                                                     self.start]
        return s

    def downhstream_igr5(self):
        pass

# url = chromosome_urls.splitlines()[0]


for url in chromosome_urls.splitlines():
    fn = bn(up(url).path)
    pth = data_dir/fn
    if pth.exists():
        local_last_mod = dt.fromtimestamp(pth.stat().st_mtime).astimezone(
            tz.utc)
    else:
        local_last_mod = dt.fromtimestamp(0, tz=tz.utc)
    response = requests.get(url, stream=True)
    remote_last_mod = parsedate_to_datetime(
        response.headers.get('last-modified') or 0)
    if local_last_mod != remote_last_mod:
        time_stamp = calendar.timegm(remote_last_mod.timetuple())
        total = int(response.headers.get('content-length'))
        with open(pth, 'wb') as f:
            for data in tqdm(response.iter_content(), total=total):
                f.write(data)
        os.utime(pth, times=(time_stamp,)*2)
    else:
        time_stamp = calendar.timegm(local_last_mod.timetuple())

    krom = SeqIO.read(pth, "gb")
    fastapath = pth.with_suffix(".fasta")
    SeqIO.write(krom, fastapath, "fasta")
    os.utime(fastapath, times=(time_stamp,)*2)
    accession = krom.id
    features = dd(list)
    for f in krom.features:
        features[f.type].append(f)
    cdslist = features["CDS"]
    genelist = []
    for cds in cdslist:
        sysname = cds.qualifiers["locus_tag"][0]
        stdname = (cds.qualifiers.get("gene") or [None]).pop()
        genelist.append(Gene(sysname,
                             stdname,
                             int(cds.location.start),
                             int(cds.location.end),
                             cds.location.strand,
                             accession,
                             fastapath))
    for i in range(1, len(genelist)-1):
        genelist[i].pred = genelist[i-1]
        genelist[i].succ = genelist[i+1]

    genelist[0].succ = genelist[1]
    genelist[-1].pred = genelist[-2]

    genes = CaseInsensitiveDict((g.sysname, g) for g in genelist)
    stdgenes = CaseInsensitiveDict()

    for key in genes:
        stdname = genes[key].stdname
        if stdname:
            stdgenes.update(((stdname, genes[key]),))

    genes.update(stdgenes)

    break



    # https://www.yeastgenome.org/locus/cdc19