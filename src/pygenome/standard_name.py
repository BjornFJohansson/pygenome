#!/usr/bin/env python
# -*- coding: utf-8 -*-
import platform as _platform

if _platform.python_version().startswith("3.8"):
    import pickle        as _pickle
else:
    import pickle5       as _pickle
import os               as _os
import re               as _re
from pygenome._pretty   import pretty_str as _ps

data_dir = _os.path.join( _os.getenv("pygenome_data_dir"), "Saccharomyces_cerevisiae")
_systematic_to_standard = _pickle.load( open(_os.path.join(data_dir, "systematic_to_standard.pickle"), "rb" ) )
_feature_list = _pickle.load( open(_os.path.join(data_dir, "feature_list.pickle"), "rb" ) )


def _standard_name(gene):

    gene = gene.upper()

    if not _re.match(r"Y[A-P](R|L)\d{3}(W|C)(-.)*", gene[:7]):
        raise ValueError("{} is not a systematic gene name.".format(gene))
    else:
        try:
            gene = _systematic_to_standard[gene]
        except KeyError:
            return None
    return _ps(gene)