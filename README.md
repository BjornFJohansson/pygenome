pygenome
========

[![image]]

[![image][1]]

[![Documentation Status]]

[![image][2]]

Pygenome provide access to the [Saccharomyces cerevisiae] genome from Python. [Genes], [promoters], [terminators], and [intergenic], 
sequences as well as the deletion [loci] created by the [genome wide deletion project] are available by their systematic names 
(like [YPR080w]) or by standard name (like [CYC1]). DNA sequences are returned as Biopython [SeqRecord] objects.

Typical usage at the [IPython] command line could look like this:

    from pygenome import sg

    sg.gene["TEF1"]
    Out[2]: yeast gene YPR080W

    sg.gene["TEF1"].cds()
    Out[3]: SeqRecord(seq=Seq('ATGGGTAAAGAGAAGTCTCACATTAACGTTGTCGTTATCGGTCATGTCGATTCT...TAA', IUPACAmbiguousDNA()), 
    id='BK006949.2', name='BK006949', description='BK006949 REGION: 700594..701970', dbxrefs=[])

    sg.gene["TEF1"].locus()
    Out[4]: SeqRecord(seq=Seq('CTTCATCGGTATCTTCGCTATATTCTTTTTAGTCGAATTTGCGGGGAGAAGATG...AAC', IUPACAmbiguousDNA()), 
    id='BK006949.2', name='BK006949', description='BK006949 REGION: 699594..702970', dbxrefs=[])

    sg.gene["TEF1"].promoter()
    Out[5]: SeqRecord(seq=Seq('ACAATGCATACTTTGTACGTTCAAAATACAATGCAGTAGATATATTTATGCATA...AAA', IUPACAmbiguousDNA()), 
    id='YPR079W_YPR080W', name='.', description='BK006949 REGION: 700015..700593', dbxrefs=[])

    sg.gene["TEF1"].deletion_locus()
    Out[6]: 'No deletion primers available!'

     sg.gene["CYC1"].deletion_locus()
    Out[7]: SeqRecord(seq=Seq('GAGGCACCAGCGTCAGCATTTTCAAAGGTGTGTTCTTCGTCAGACATGTTTTAG...GTG', IUPACAmbiguousDNA()), 
    id='yjr048w::KanMX4 locus with 1000 bp up and 1000 bp downstream DNA', name='yjr048w::KanMX4', description='<unknown description>', dbxrefs=[])

|         |      |         |
|---------|------|---------|
| version | date | comment |

  [image]: https://travis-ci.org/BjornFJohansson/pygenome.svg
  [![image]]: https://travis-ci.org/BjornFJohansson/pygenome
  [1]: https://coveralls.io/repos/BjornFJohansson/pygenome/badge.svg?branch=master
  [![image][1]]: https://coveralls.io/r/BjornFJohansson/pygenome?branch=master
  [Documentation Status]: https://readthedocs.org/projects/pygenome/badge/?version=latest
  [![Documentation Status]]: https://readthedocs.org/projects/pygenome/?badge=latest
  [2]: https://badge.fury.io/py/pygenome.svg
  [![image][2]]: http://badge.fury.io/py/pygenome
  [Saccharomyces cerevisiae]: https://microbewiki.kenyon.edu/index.php/Saccharomyces_cerevisiae
  [Genes]: http://en.wikipedia.org/wiki/Gene
  [promoters]: http://en.wikipedia.org/wiki/Promoter_(genetics)
  [terminators]: http://en.wikipedia.org/wiki/Terminator_(genetics)
  [intergenic]: http://en.wikipedia.org/wiki/Intergenic_region
  [loci]: http://en.wikipedia.org/wiki/Locus_(genetics)
  [genome wide deletion project]: http://www-sequence.stanford.edu/group/yeast_deletion_project/deletions3.html
  [YPR080w]: http://www.yeastgenome.org/locus/S000006284/overview
  [CYC1]: http://www.yeastgenome.org/locus/S000003809/overview
  [SeqRecord]: http://biopython.org/wiki/SeqRecord
  [IPython]: http://ipython.org/
