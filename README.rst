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

.. image:: https://pypip.in/download/pygenome/badge.svg
    :target: https://pypi.python.org/pypi/pygenome/
    :alt: Downloads
    
.. image:: https://pypip.in/version/pygenome/badge.svg
    :target: https://pypi.python.org/pygenome/pydna/
    :alt: Latest Version

.. image:: https://pypip.in/wheel/pydna/badge.svg
    :target: https://pypi.python.org/pypi/pydna/
    :alt: Wheel Status

Pygenome provide a module for accessing the Saccharomyces cerevisiae genome from 
Python. Genes, promoters, terminators and intergenic
sequences are available by systematic names (like YPR080w) or by standard name.
DNA sequences are returned as Biopython SeqRecord objects. Biopython, percache and
appdirs are required for installation.

Typical usage at the command line could look like this::

    >>> from pygenome import sg            
    >>> sg.locus("TEF1")
    SeqRecord(seq=Seq('CTTCATCGGTATCTTCGCTATATTCTTTTTAGTCGAATTTGCGGGGAGAAGATG...AAC', IUPACAmbiguousDNA()), id='BK006949.2', name='BK006949', description='TPA: Saccharomyces cerevisiae S288c chromosome XVI.', dbxrefs=[])
    >>> sg.locus("YPR080w")
    SeqRecord(seq=Seq('CTTCATCGGTATCTTCGCTATATTCTTTTTAGTCGAATTTGCGGGGAGAAGATG...AAC', IUPACAmbiguousDNA()), id='BK006949.2', name='BK006949', description='TPA: Saccharomyces cerevisiae S288c chromosome XVI.', dbxrefs=[])
    >>> a=sg.locus("TEF1")
    >>> b=sg.locus("YPR080w")
    >>> str(a.seq) == str(b.seq)
    True
    >>> prom=sg.promoter("TEF1")
    >>> prom
    SeqRecord(seq=Seq('ACAATGCATACTTTGTACGTTCAAAATACAATGCAGTAGATATATTTATGCATA...AAA', IUPACAmbiguousDNA()), id='BK006949.2', name='BK006949', description='TPA: Saccharomyces cerevisiae S288c chromosome XVI.', dbxrefs=[])
    >>> term=sg.terminator("TEF1")
    >>> term
    SeqRecord(seq=Seq('GGAGATTGATAAGACTTTTCTAGTTGCATATCTTTTATATTTAAATCTTATCTA...CAG', IUPACAmbiguousDNA()), id='BK006949.2', name='BK006949', description='TPA: Saccharomyces cerevisiae S288c chromosome XVI.', dbxrefs=[])
    
    
NEWS
====


=======   ========== =============================================================
version   date       comment
=======   ========== =============================================================
0.5.0     2015-03-03 Documentation, automatic build,test and deployment

0.0.6     2014-06-17 Bugfix

0.0.5     2014-06-14 Simpler api (see example above)

0.0.1     2013-08-01 first release
=======   ========== =============================================================


System Requirements
===================

- `Python 2.7 <http://www.python.org>`_.

- `biopython <http://pypi.python.org/pypi/biopython>`_.

- `percache  <http://pypi.python.org/pypi/percache>`_.

- `appdirs <http://pypi.python.org/pypi/appdirs>`_.

Python 2.x
----------

Versions other than 2.7 has not been tried with this software.
Version 2.7.3 was used to build the distribution.

Python 3.x
----------

This code has not been tested with python 3.

Installation
============

Source
------

Open the pydna source code directory (containing the setup.py file) in
terminal and type:

    sudo python setup.py install <enter>

If you need to do additional configuration, e.g. changing the base
directory, please type `python setup.py`, or see the documentation for
Setuptools.


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
