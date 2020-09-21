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
_standard_to_systematic = _pickle.load( open(_os.path.join(data_dir, "standard_to_systematic.pickle"), "rb" ) )
_feature_list = _pickle.load( open(_os.path.join(data_dir, "feature_list.pickle"), "rb" ) )


def _systematic_name(gene):

    gene = gene.upper()

    if _re.match(r"Y[A-P](R|L)\d{3}(W|C)(-.)*", gene[:7]) and gene in _feature_list:
        return _ps(gene)
    else:
        try:
            gene = _standard_to_systematic[gene]
        except KeyError:
            raise KeyError("gene {} does not exist".format(gene))
    return _ps(gene)
