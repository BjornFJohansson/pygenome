#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on %(date)s

@author: %(username)s
"""

from pyfaidx import Fasta
chromosome = Fasta('chr01.fsa')

chromosome.keys()

chromosome['tpg|BK006935.2|'][1807:2169]



