pygenome
========

![image](https://travis-ci.org/BjornFJohansson/pygenome.svg%20%0A%20:target:%20https://travis-ci.org/BjornFJohansson/pygenome)

![image](https://coveralls.io/repos/BjornFJohansson/pygenome/badge.svg?branch=master%20%0A%20:target:%20https://coveralls.io/r/BjornFJohansson/pygenome?branch=master)

![image](https://readthedocs.org/projects/pygenome/badge/?version=latest%0A%20:target:%20https://readthedocs.org/projects/pygenome/?badge=latest%0A%20:alt:%20Documentation%20Status)

![image](https://badge.fury.io/py/pygenome.svg%0A%20:target:%20http://badge.fury.io/py/pygenome)

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

| ver   | date       | comment                                             |
|-------|------------|-----------------------------------------------------|
| 0.9.5 | 2017-01-01 | Python 3 release                                    |
| 0.9.0 | 2015-05-01 | Changed interface to a more object oriented style   |
| 0.5.0 | 2015-03-03 | Documentation, automatic build, test and deployment |
| 0.0.6 | 2014-06-17 | Bugfix                                              |
| 0.0.5 | 2014-06-14 | Simpler api (see example above)                     |
| 0.0.1 | 2013-08-01 | first release                                       |

System Requirements
-------------------

-   [Python 3](http://www.python.org) (0.9.0 was the last to support Python 2.7.)
-   [biopython](http://pypi.python.org/pypi/biopython)
-   [percache](http://pypi.python.org/pypi/percache)
-   [appdirs](http://pypi.python.org/pypi/appdirs)

Installation
------------



The second best way of installing pygenome is by using
[pip](https://packaging.python.org/en/latest/installing.html#installing-from-pypi)

> sudo pip install pygenome

### Source Code Repository

pydna source code is hosted on [Github](https://github.com/BjornFJohansson/pygenome)
