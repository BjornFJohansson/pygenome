#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


"""

# http://sgd-archive.yeastgenome.org/sequence/S288C_reference/NCBI_genome_source

import os
from email.utils import parsedate_to_datetime
from tqdm import tqdm
import requests
from pathlib import Path
from urllib.parse import urlparse
from os.path import basename
from datetime import timedelta

data_path = Path("data")

# urls for each of the sixteen chromosomes
file_path = data_path / Path("chromosome_urls.txt")

# Read the file and store each line as an element in a list
chromosome_urls = file_path.read_text(encoding="utf-8").splitlines()


for url in chromosome_urls:

    r = requests.get(url, stream=True)

    remote_last_changed = parsedate_to_datetime(r.headers.get('last-modified') or 0)

    gb_pth = data_path / Path(basename(urlparse(url).path))

    if not gb_pth.exists() or gb_pth.stat().st_mtime < remote_last_changed.timestamp():
        total = int(r.headers.get('content-length'))
        with open(gb_pth, 'wb') as f:
            for data in tqdm(r.iter_content(),
                             desc=str(gb_pth),
                             total=total):
                f.write(data)
        os.utime(gb_pth, times=(remote_last_changed.timestamp(),
                                remote_last_changed.timestamp()))
    else:
        print(f"{gb_pth.name} is up to date.")

#     if not fa_pth.exists() or gb_pth.stat().st_mtime > fa_pth.stat().st_mtime:
#         krom = SeqIO.read(gb_pth, "gb")

#         pkl = pickle.dumps(tuple(krom.features), 5)
#         pkl = _pickletools.optimize(pkl)




# nd = []
# pth = _Path("ORFs_not_available.txt")
# with open(pth, 'rt') as csvfile:
#     rd = csv.reader(csvfile, delimiter='\t')
#     next(rd)
#     next(rd)
#     next(rd)
#     next(rd)
#     next(rd)
#     not_done = {}
#     for line_ in rd:
#         nd.append([x.strip() for x in line_][0])

# pkl = pickle.dumps(tuple(nd), 5)
# pkl = _pickletools.optimize(pkl)
# zipObj.writestr(str(pth.with_suffix(".pickle")),
#                 pkl,
#                 compress_type=ZIP_DEFLATED,
#                 compresslevel=9)

# # Deletion_primers_PCR_sizes_corrected
# # Some ups45 anneal to the plasmid but are wrong.
# # changed the sequences below
# # GATGTCCACGAGCTCTCT -> GATGTCCACGAGGTCTCT

# pth = _Path("Deletion_primers_PCR_sizes_corrected.txt")
# primers = dd(list)

# with open(pth, 'rt') as csvfile:
#     rd = csv.reader(csvfile, delimiter='\t')
#     next(rd)
#     next(rd)
#     for line_ in rd:
#         v = tuple([x.strip() for x in line_])
#         primers[v[1]].append(v)


# primers = {k: tuple(v) for k, v in primers.items()}
# pkl = pickle.dumps(primers, 5)
# pkl = _pickletools.optimize(pkl)
# zipObj.writestr("Deletion_primers_PCR_sizes.pickle",
#                 pkl,
#                 compress_type=ZIP_DEFLATED,
#                 compresslevel=9)


# # 0	rec_num
# # 1	ORF_name
# # 2	deletion_alias
# # 3	essential
# # 4	A_confirmation_primer_sequence
# # 5	B_confirmation_primer_sequence
# # 6	C_confirmation_primer_sequence
# # 7	D_confirmation_primer_sequence
# # 8	UPTAG_primer_sequence
# # 9	DNTAG_primer_sequence
# # 10	UPstream45_primer_sequence
# # 11	DNstream45_primer_sequence
# # 12	UPstream90_primer_sequence
# # 13	DNstream90_primer_sequence
# # 14	AB_wt_PCR
# # 15	AkanB_del_PCR
# # 16	CD_wt_PCR
# # 17	DkanC_del_PCR
# # 18	AD_wt
# # 19	AD_del
# # 20	AB_del_PCR
# # 21	CD_del_PCR
# # 22	UPTAG_sequence_20mer
# # 23	DNTAG_sequence_20mer

# pth = _Path("yeastGFPOligoSequence.txt")
# primers = dd(list)

# with open(pth, 'rt') as csvfile:
#     rd = csv.reader(csvfile, delimiter='\t')
#     for line_ in rd:
#         v = tuple([x.strip() for x in line_])
#         primers[v[0]] = v

# pkl = pickle.dumps(primers, 5)
# pkl = _pickletools.optimize(pkl)
# zipObj.writestr("yeastGFPOligoSequence.pickle",
#                 pkl,
#                 compress_type=ZIP_DEFLATED,
#                 compresslevel=9)




# # close the Zip File
# zipObj.close()
