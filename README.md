# ![icon](https://raw.githubusercontent.com/bjornFJohansson/pygenome/master/docs/SuperYeast.jpg) pygenome

[![codecov](https://codecov.io/gh/BjornFJohansson/pygenome/branch/master/graph/badge.svg)](https://codecov.io/gh/BjornFJohansson/pygenome)
[![Documentation Status](https://readthedocs.org/projects/pygenome/badge/?version=latest)](http://pygenome.readthedocs.io/en/latest/?badge=latest)
[![PyPI version](https://badge.fury.io/py/pygenome.svg)](https://badge.fury.io/py/pygenome)
[![Anaconda-Server Badge](https://anaconda.org/bjornfjohansson/pygenome/badges/version.svg)](https://anaconda.org/bjornfjohansson/pygenome)

Harness the awesome power of yeast genetics through python! Pygenome provide access to the [Saccharomyces cerevisiae](https://microbewiki.kenyon.edu/index.php/Saccharomyces_cerevisiae)
genome from Python. [Genes](http://en.wikipedia.org/wiki/Gene),
[promoters](http://en.wikipedia.org/wiki/Promoter_(genetics)),
[terminators](http://en.wikipedia.org/wiki/Terminator_(genetics)), and
[intergenic](http://en.wikipedia.org/wiki/Intergenic_region), sequences
as well as the deletion [loci](http://en.wikipedia.org/wiki/Locus_(genetics)) created by the
[genome wide deletion project](http://www-sequence.stanford.edu/group/yeast_deletion_project/deletions3.html)
are available by their systematic names (like [YPR080w](http://www.yeastgenome.org/locus/S000006284/overview)) or by
standard name (like [CYC1](http://www.yeastgenome.org/locus/S000003809/overview)). DNA
sequences are returned as Biopython
[SeqRecord](http://biopython.org/wiki/SeqRecord) objects. Thanks to [SGD](http://www.yeastgenome.org/) for letting me use the SuperYeast logotype above.

Typical usage at the [IPython](http://ipython.org/) command line could look like this:

```python
from pygenome import sg
```


```python
mygene = sg.stdgene["XKS1"]
```


```python
mygene
```




<a href='http://www.yeastgenome.org/locus/YGR194C' target='_blank'>Gene XKS1/YGR194C</a>




```python
mygene.short_description
```




    Xylulokinase; converts D-xylulose and ATP to xylulose 5-phosphate and ADP; rate limiting step in fermentation of xylulose; required for xylose fermentation by recombinant S. cerevisiae strains




```python
sg.sysgene["YGR194C"]
```




<a href='http://www.yeastgenome.org/locus/YGR194C' target='_blank'>Gene XKS1/YGR194C</a>




```python
mygene.cds
```




    SeqRecord(seq=Seq('ATGTTGTGTTCAGTAATTCAGAGACAGACAAGAGAGGTTTCCAACACAATGTCT...TAA'), id='BK006941.2', name='BK006941', description='BK006941 REGION: complement(887876..886072)', dbxrefs=[])




```python
mygene.locus()
```




    SeqRecord(seq=Seq('ATCCTGCTGTAGTTATGGCACTAAAGTTTTTTTGTAAATCTTTTTATATGTTAA...GAA'), id='BK006941.2', name='BK006941', description='BK006941 REGION: complement(888876..885072)', dbxrefs=[])




```python
mygene.promoter
```




    SeqRecord(seq=Seq('ATGATGATCCTGCTGTAGTTATGGCACTAAAGTTTTTTTGTAAATCTTTTTATA...TTA'), id='YGR195W_YGR194C', name='.', description='BK006941.2 REGION: complement(887876..888881)', dbxrefs=[])




```python
mygene.terminator
```




    SeqRecord(seq=Seq('AATATGTTTGAATAATTTATCATGCCCTGACAAGTACACACAAACACAGACACA...AAA'), id='YGR194C_YGR195W', name='.', description='Intergenic sequence between upstream gene YGR194C and downstream gene Gene PDX1/YGR193C', dbxrefs=[])




```python
mygene.downstream_gene
```




<a href='http://www.yeastgenome.org/locus/YGR193C' target='_blank'>Gene PDX1/YGR193C</a>




```python
mygene.upstream_gene
```




<a href='http://www.yeastgenome.org/locus/YGR195W' target='_blank'>Gene SKI6/YGR195W</a>




```python
mygene.deletion_locus
```




    SeqRecord(seq=Seq('ATCCTGCTGTAGTTATGGCACTAAAGTTTTTTTGTAAATCTTTTTATATGTTAA...GAA'), id='ygr194c::KanMX4 locus with 1000 bp up and 1000 bp downstream DNA', name='ygr194c::KanMX4', description='description', dbxrefs=[])



http://www-sequence.stanford.edu/group/yeast_deletion_project/downloads.html

| ver   | date       | comment                                             |
|-------|------------|-----------------------------------------------------|
| 2.0.0 | 2017-09-02 | split sg.gene dict into sg.stdgene and sg.sysgene   |
| 1.0.0 | 2017-03-24 | Internal stuff, automativ build & test              |
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

- [Python 3.6 or 3.7](http://www.python.org) (pygenome version 0.9.0 was the last to support Python 2.7.)
- [pydna](http://pypi.python.org/pypi/pydna)
- [requests](http://pypi.python.org/pypi/requests)
- [appdirs](https://pypi.python.org/pypi/appdirs)

## Install with pip

The second best way of installing pygenome is by using
[pip](https://packaging.python.org/en/latest/installing.html#installing-from-pypi)

    sudo pip install pygenome

### Source Code Repository

pydna source code is hosted on [Github](https://github.com/BjornFJohansson/pygenome).
