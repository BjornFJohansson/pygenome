#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright 2013, 2014, 2015 by Björn Johansson.  All rights reserved.
# This code is part of the pygenome distribution and governed by its
# license.  Please see the LICENSE.txt file that should have been included
# as part of this package.


__version__      = "0.06"
__date__         = "2014-06-17"
__author__       = "Björn Johansson"
__copyright__    = "Copyright 2013,2014, Björn Johansson"
__credits__      = ["Björn Johansson"]
__license__      = "BSD"
__maintainer__   = "Björn Johansson"
__email__        = "bjorn_johansson@bio.uminho.pt"
__status__       = "Development" # "Production" #"Prototype"

from saccharomyces_cerevisiae import saccharomyces_cerevisiae_genome

sg = saccharomyces_cerevisiae_genome()

#del saccharomyces_cerevisiae

del saccharomyces_cerevisiae_genome

