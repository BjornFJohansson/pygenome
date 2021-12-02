#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Access to the Saccharomyces cerevisiae genome from Python.

Sequences can be accessed as Bio.SeqRecord objects provided by Biopython.
"""
import os
import calendar
import datetime
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

from Bio import SeqIO
from pydna._pretty import pretty_str as ps
from requests.structures import CaseInsensitiveDict

if sys.version_info >= (3, 8):
    import pickle
else:
    import pickle5 as pickle

data_dir = Path(os.getenv("pygenome_data_dir"))/"Saccharomyces_cerevisiae"


# set the ch_urls and ch_file_names variables
with open(data_dir/"chromosome_urls.py", "r") as f:
    exec(f.read(), globals(), locals())
assert chromosome_urls

ac_dct = {1: "{}..{}",
          -1: "complement({}..{})"}

# url = chromosome_urls.splitlines()[0]

for url in chromosome_urls.splitlines():
    fn = bn(up(url).path)
    pth = data_dir/fn
    if pth.exists():     
        local_last_mod = datetime.datetime.fromtimestamp(pth.stat().st_mtime).astimezone(datetime.timezone.utc)
    else:
        local_last_mod = datetime.datetime.fromtimestamp(0, tz=datetime.timezone.utc)
    response = requests.get(url, stream=True)
    remote_last_mod = parsedate_to_datetime(response.headers.get('last-modified') or 0)
    if local_last_mod != remote_last_mod:
        remote_last_mod_time_stamp = calendar.timegm(remote_last_mod.timetuple())
        total = int(response.headers.get('content-length'))
        with open(pth, 'wb') as f:
            for data in tqdm(response.iter_content(), total=total):
                f.write(data)
        os.utime(pth, times=(remote_last_mod_time_stamp,)*2)
    krom = SeqIO.read(pth, "gb")
    fastapath = pth.with_suffix(".fasta")
    SeqIO.write(krom, fastapath, "fasta")
    features = dd(list)
    for f in krom.features:
        features[f.type].append(f)

    # del response
    # del krom
    # CDS = features["CDS"][1]    

    from itertools import zip_longest
    
    for PRD, CDS, SUC in zip_longest(*[iter(features["CDS"])] * 3, fillvalue=None):
        pass
        sysname = CDS.qualifiers["locus_tag"][0]
        stdname = CDS.qualifiers.get("gene") or None
        
        
        
        try:
            description = f.qualifiers["note"][0]
        except KeyError:
            description = f.qualifiers["product"][0]
        descriptions.append(description)








