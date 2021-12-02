#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""docstring."""

import time
import os
import logging
from logging import handlers
from pathlib import Path
import shutil
import zipfile
from pkg_resources import resource_filename

# create logger
loglevel = os.environ["pygenome_loglevel"]
logdir = os.environ["pygenome_log_dir"]
logger = logging.getLogger("pygenome")
logger.setLevel(loglevel)
hdlr = handlers.RotatingFileHandler(os.path.join(logdir, 'pygenome.log'),
                                    mode='a',
                                    maxBytes=10*1024*1024,
                                    backupCount=10,
                                    encoding='utf-8')
formatter = logging.Formatter(('{asctime} {levelname}'
                               ' {funcName} {message}'), style='{')
hdlr.setFormatter(formatter)
logger.addHandler(hdlr)
logger.info('Logger started.')

data_dir = Path(os.environ["pygenome_data_dir"])/"Saccharomyces_cerevisiae"
data_dir.mkdir(parents=True, exist_ok=True)
logger.info("data_dir set to:  %s", data_dir)

datafile = Path(resource_filename("pygenome",
                ("saccharomyces_cerevisiae/"
                 "Saccharomyces_cerevisiae.zip")))


with zipfile.ZipFile(datafile, "r") as zf:
    for zi in zf.infolist():
        pth = data_dir/zi.filename
        if not pth.exists():
            zf.extract(zi, path=data_dir)
            date_time = time.mktime(zi.date_time + (0, 0, -1))
            os.utime(pth, (date_time, date_time))


config_dir = Path(os.environ["pygenome_config_dir"])

sf = config_dir/"settings_Saccharomyces_cerevisiae.py"

if not sf.exists():
    shutil.copy(data_dir/"settings_Saccharomyces_cerevisiae.py", config_dir)


pickles = """\
primers.pickle
not_done.pickle
feature_list.pickle
standard_to_systematic.pickle
systematic_to_standard.pickle
systematic_to_genbank_accession.pickle
systematic_to_description.pickle
""".splitlines()


if not all([os.path.exists(data_dir/pickle) for pickle in pickles]):
    logger.info("pickles start.")
    from pygenome.saccharomyces_cerevisiae._pickle_primers import pickle_primers
    from pygenome.saccharomyces_cerevisiae._pickle_primers import pickle_orfs_not_deleted
    pickle_primers()
    pickle_orfs_not_deleted()
    from pygenome._pickle_lists import _pickle_lists
    from pygenome._pickle_genes import _pickle_genes
    pickle_lists()
    pickle_genes()
    logger.info("pickles done.")
