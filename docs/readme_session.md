```python
from pygenome import saccharomyces_cerevisiae as sg
```


```python
mygene = sg.stdgenes["XKS1"]
```


```python
mygene
```




<a href='http://www.yeastgenome.org/locus/YGR194C' target='_blank'>Gene XKS1/YGR194C</a>




```python
mygene.short_description()
```




    Xylulokinase; converts D-xylulose and ATP to xylulose 5-phosphate and ADP; rate limiting step in fermentation of xylulose; required for xylose fermentation by recombinant S. cerevisiae strains




```python
sg.sysgenes["YGR194C"]
```




<a href='http://www.yeastgenome.org/locus/YGR194C' target='_blank'>Gene XKS1/YGR194C</a>




```python
mygene.cds()
```




    Dseqrecord(-1803)




```python
mygene.locus()
```




    Dseqrecord(-3803)




```python
mygene.promoter()
```




    Dseqrecord(-1006)




```python
mygene.terminator()
```




    Dseqrecord(-331)




```python
mygene.downstream_gene()
```




<a href='http://www.yeastgenome.org/locus/YGR193C' target='_blank'>Gene PDX1/YGR193C</a>




```python
mygene.upstream_gene()
```




<a href='http://www.yeastgenome.org/locus/YGR195W' target='_blank'>Gene SKI6/YGR195W</a>




```python
mygene.deletion_cassettes()
```




    [Cassette(outer_cassette=Amplicon(1671), inner_cassette=Amplicon(1617), UPTAG=UPTAG_primer_YGR194C 74-mer:5'-GAGATTAGTACTTTA..GAC-3', DNTAG=DNTAG_primer_YGR194C 74-mer:5'-TTATTCAAACATATT..TCG-3', UPstream45=UPstream45_YGR194C 45-mer:5'-CCCTCTCGAGAAAAA..ATG-3', DNstream45=DNstream45_YGR194C 45-mer:5'-TGTGTGTACTTGTCA..TTA-3')]




```python
mygene.deletion_loci()
```




    [Contig(-3587)]


