========
pygenome
========

.. image:: https://travis-ci.org/BjornFJohansson/pygenome.svg 
    :target: https://travis-ci.org/BjornFJohansson/pygenome
    
.. image:: https://coveralls.io/repos/BjornFJohansson/pygenome/badge.svg?branch=master 
    :target: https://coveralls.io/r/BjornFJohansson/pygenome?branch=master
  
.. image:: https://readthedocs.org/projects/pygenome/badge/?version=latest
    :target: https://readthedocs.org/projects/pygenome/?badge=latest
    :alt: Documentation Status
    
.. image:: https://badge.fury.io/py/pygenome.svg
    :target: http://badge.fury.io/py/pygenome

Pygenome provide access to the `Saccharomyces cerevisiae <https://microbewiki.kenyon.edu/index.php/Saccharomyces_cerevisiae>`_ genome from 
Python. `Genes <http://en.wikipedia.org/wiki/Gene>`_, `promoters <http://en.wikipedia.org/wiki/Promoter_(genetics)>`_,
`terminators <http://en.wikipedia.org/wiki/Terminator_(genetics)>`_, and 
`intergenic <http://en.wikipedia.org/wiki/Intergenic_region>`_,
sequences as well as the deletion `loci <http://en.wikipedia.org/wiki/Locus_(genetics)>`_ created by the `genome wide deletion project <http://www-sequence.stanford.edu/group/yeast_deletion_project/deletions3.html>`_ 
are available by their systematic names (like `YPR080w <http://www.yeastgenome.org/locus/S000006284/overview>`_) or by 
standard name (like `CYC1 <http://www.yeastgenome.org/locus/S000003809/overview>`_).
DNA sequences are returned as Biopython `SeqRecord <http://biopython.org/wiki/SeqRecord>`_ objects.

Typical usage at the `IPython <http://ipython.org/>`_ command line could look like this::


    from pygenome import sg

    sg.gene["TEF1"]
    Out[2]: yeast gene YPR080W

    sg.gene["TEF1"].cds()
    Out[3]: SeqRecord(seq=Seq('ATGGGTAAAGAGAAGTCTCACATTAACGTTGTCGTTATCGGTCATGTCGATTCT...TAA', IUPACAmbiguousDNA()), id='BK006949.2', name='BK006949', description='BK006949 REGION: 700594..701970', dbxrefs=[])

    sg.gene["TEF1"].locus()
    Out[4]: SeqRecord(seq=Seq('CTTCATCGGTATCTTCGCTATATTCTTTTTAGTCGAATTTGCGGGGAGAAGATG...AAC', IUPACAmbiguousDNA()), id='BK006949.2', name='BK006949', description='BK006949 REGION: 699594..702970', dbxrefs=[])

    sg.gene["TEF1"].promoter()
    Out[5]: SeqRecord(seq=Seq('ACAATGCATACTTTGTACGTTCAAAATACAATGCAGTAGATATATTTATGCATA...AAA', IUPACAmbiguousDNA()), id='YPR079W_YPR080W', name='.', description='BK006949 REGION: 700015..700593', dbxrefs=[])

    sg.gene["TEF1"].deletion_locus()
    Out[6]: 'No deletion primers available!'

     sg.gene["CYC1"].deletion_locus()
    Out[7]: SeqRecord(seq=Seq('GAGGCACCAGCGTCAGCATTTTCAAAGGTGTGTTCTTCGTCAGACATGTTTTAG...GTG', IUPACAmbiguousDNA()), id='yjr048w::KanMX4 locus with 1000 bp up and 1000 bp downstream DNA', name='yjr048w::KanMX4', description='<unknown description>', dbxrefs=[])




=======   ========== =============================================================
version   date       comment
=======   ========== =============================================================
0.9.0     2015-05-01 Changed interface to a more object oriented style

0.5.0     2015-03-03 Documentation, automatic build, test and deployment

0.0.6     2014-06-17 Bugfix

0.0.5     2014-06-14 Simpler api (see example above)

0.0.1     2013-08-01 first release
=======   ========== =============================================================


System Requirements
===================

- `Python 3.5 <http://www.python.org>`_

- `biopython <http://pypi.python.org/pypi/biopython>`_

- `percache  <http://pypi.python.org/pypi/percache>`_

- `appdirs <http://pypi.python.org/pypi/appdirs>`_

Python 2.7
----------

Version 0.9.0 was the last one to support Python 2.7.

Installation
============

Source
------

The best way of installing pygenome is by using `pip <https://packaging.python.org/en/latest/installing.html#installing-from-pypi>`_

    sudo pip install pygenome


Source Code Repository
----------------------

pydna source code is hosted on Github:

https://github.com/BjornFJohansson/pygenome


Distribution Structure
======================

README.rst          -- This file.

LICENSE.txt         -- What you can do with the code.

MANIFEST.in         -- Tells distutils what files to distribute

setup.py            -- Installation file.

pygenome/           -- The actual code.

docs/               -- Documentation.

scripts/            -- Miscellaneous perhaps useful scripts

tests/              -- Testing code
