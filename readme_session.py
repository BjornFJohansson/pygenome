#!/usr/bin/env python
# -*- coding: utf-8 -*-
from pygenome import sg
mygene = sg.stdgene["XKS1"]
mygene
mygene.short_description
sg.sysgene["YGR194C"]
mygene.cds
mygene.locus()
mygene.promoter
mygene.terminator
mygene.downstream_gene
mygene.upstream_gene
mygene.deletion_locus


