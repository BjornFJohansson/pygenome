# pygenome

[![Build Status](https://travis-ci.org/BjornFJohansson/pygenome.svg?branch=master)](https://travis-ci.org/BjornFJohansson/pygenome)[![Coverage Status](https://coveralls.io/repos/github/BjornFJohansson/pygenome/badge.svg?branch=master)](https://coveralls.io/github/BjornFJohansson/pygenome?branch=master)[![Documentation Status](https://readthedocs.org/projects/pygenome/badge/?version=latest)](http://pygenome.readthedocs.io/en/latest/?badge=latest)
                

Pygenome provide access to the [Saccharomyces cerevisiae](https://microbewiki.kenyon.edu/index.php/Saccharomyces_cerevisiae)
genome from Python. [Genes](http://en.wikipedia.org/wiki/Gene),
[promoters](http://en.wikipedia.org/wiki/Promoter_(genetics)),
[terminators](http://en.wikipedia.org/wiki/Terminator_(genetics)), and
[intergenic](http://en.wikipedia.org/wiki/Intergenic_region), sequences
as well as the deletion [loci](http://en.wikipedia.org/wiki/Locus_(genetics)) created by the
[genome wide deletion project](http://www-sequence.stanford.edu/group/yeast_deletion_project/deletions3.html)
are available by their systematic names (like [YPR080w](http://www.yeastgenome.org/locus/S000006284/overview)) or by
standard name (like [CYC1](http://www.yeastgenome.org/locus/S000003809/overview)). DNA
sequences are returned as Biopython
[SeqRecord](http://biopython.org/wiki/SeqRecord) objects.

Typical usage at the [IPython](http://ipython.org/) command line could look like this:

    from pygenome import sg

    sg.gene["TEF1"]
    Out[2]: yeast gene YPR080W

    sg.gene["TEF1"].cds
    Out[3]: SeqRecord(seq=Seq('ATGGGTAAAGAGAAGTCTCACATTAACGTTGTCGTTATCGGTCATGTCGATTCT...TAA',
            IUPACAmbiguousDNA()), id='BK006949.2', name='BK006949', description='BK006949 REGION: 700594..701970', dbxrefs=[])

    sg.gene["TEF1"].locus()
    Out[4]: SeqRecord(seq=Seq('CTTCATCGGTATCTTCGCTATATTCTTTTTAGTCGAATTTGCGGGGAGAAGATG...AAC',
            IUPACAmbiguousDNA()), id='BK006949.2', name='BK006949', description='BK006949 REGION: 699594..702970', dbxrefs=[])

    sg.gene["TEF1"].promoter
    Out[5]: SeqRecord(seq=Seq('ACAATGCATACTTTGTACGTTCAAAATACAATGCAGTAGATATATTTATGCATA...AAA',
            IUPACAmbiguousDNA()), id='YPR079W_YPR080W', name='.', description='BK006949 REGION: 700015..700593', dbxrefs=[])

    sg.gene["TEF1"].deletion_locus
    Out[6]: 'No deletion primers available!'

     
    Out[7]: SeqRecord(seq=Seq('GAGGCACCAGCGTCAGCATTTTCAAAGGTGTGTTCTTCGTCAGACATGTTTTAG...GTG',
            IUPACAmbiguousDNA()), id='yjr048w::KanMX4 locus with 1000 bp up and 1000 bp downstream DNA',
            name='yjr048w::KanMX4', description='<unknown description>', dbxrefs=[])

| ver   | date       | comment                                             |
|-------|------------|-----------------------------------------------------|
| 0.9.5 | 2017-01-01 | Python 3 release                                    |
| 0.9.0 | 2015-05-01 | Changed interface to a more object oriented style   |
| 0.5.0 | 2015-03-03 | Documentation, automatic build, test and deployment |
| 0.0.6 | 2014-06-17 | Bugfix                                              |
| 0.0.5 | 2014-06-14 | Simpler api (see example above)                     |
| 0.0.1 | 2013-08-01 | first release                                       |

## Installation using conda on Anaconda

The absolutely best way of installing and using pygenome is to use the 
free [Anaconda](https://store.continuum.io/cshop/anaconda) or [Miniconda](http://conda.pydata.org/miniconda.html) python distributions.

Anaconda is a large download (about 400 Mb) while Miniconda is about 40-50 Mb. 

Once Anaconda (or Miniconda) is installed, the conda package manager can be used to install pygenome 
from the [BjornFJohansson](https://anaconda.org/bjornfjohansson) package channel.

The first step is to add the channel by typing the command below followed by return:

    conda config --append channels BjornFJohansson

Then pygenome can be installed by typing the command below followed by return:

    conda install pygenome

This works on Windows, MacOSX and Linux, and installs all necessary dependencies automatically.

## Requirements

- [Python 3](http://www.python.org) (0.9.0 was the last to support Python 2.7.)
- [pydna](http://pypi.python.org/pypi/pydna)
- [percache](http://pypi.python.org/pypi/percache)
- [appdirs]()

## Install with pip

The second best way of installing pygenome is by using
[pip](https://packaging.python.org/en/latest/installing.html#installing-from-pypi)

    sudo pip install pygenome

### Source Code Repository

pydna source code is hosted on [Github](https://github.com/BjornFJohansson/pygenome).
