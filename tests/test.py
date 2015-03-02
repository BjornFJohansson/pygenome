#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
test pygenome
'''
import nose
from pygenome import sg as sc

def test_TEF1():
    tp = 'ACAATGCATACTTTGTACGTTCAAAATACAATGCAGTAGATATATTTATGCATATTACATATAATACATATCACATAGGAAGCAACAGGCGCGTTGGACTTTTAATTTTCGAGGACCGCGAATCCTTACATCACACCCAATCCCCCACAAGTGATCCCCCACACACCATAGCTTCAAAATGTTTCTACTCCTTTTTTACTCTTCCAGATTTTCTCGGACTCCGCGCATCGCCGTACCACTTCAAAACACCCAAGCACAGCATACTAAATTTCCCCTCTTTCTTCCTCTAGGGTGTCGTTAATTACCCGTACTAAAGGTTTGGAAAAGAAAAAAGAGACCGCCTCGTTTCTTTTTCTTCGTCGAAAAAGGCAATAAAAATTTTTATCACGTTTCTTTTTCTTGAAAATTTTTTTTTTTGATTTTTTTCTCTTTCGATGACCTCCCATTGATATTTAAGTTAATAAACGGTCTTCAATTTCTCAAGTTTCAGTTTCATTTTTCTTGTTCTATTACAACTTTTTTTACTTCTTGCTCATTAGAAAGAAAGCATAGCAATCTAATCTAAGTTTTAATTACAAA'
    s1 = sc.intergenic_sequence("YPR079W", "TEF1")
    s2 = sc.promoter("TEF1")
    assert tp.lower() == str(s1.seq).lower() == str(s2.seq).lower()

def test_TEF1_genbank_codes():
    assert  sc.intergenic_sequence_genbank_accession("YPR079W", "TEF1") == 'BK006949.2 REGION: 700015..700593'
    assert  sc.promoter_genbank_accession("TEF1") == 'BK006949.2 REGION: 700015..700593'

def test_TPI1():
    tp = 'TGTTTAAAGATTACGGATATTTAACTTACTTAGAATAATGCCATTTTTTTGAGTTATAATAATCCTACGTTAGTGTGAGCGGGATTTAAACTGTGAGGACCTTAATACATTCAGACACTTCTGCGGTATCACCCTACTTATTCCCTTCGAGATTATATCTAGGAACCCATCAGGTTGGTGGAAGATTACCCGTTCTAAGACTTTTCAGCTTCCTCTATTGATGTTACACCTGGACACCCCTTTTCTGGCATCCAGTTTTTAATCTTCAGTGGCATGTGAGATTCTCCGAAATTAATTAAAGCAATCACACAATTCTCTCGGATACCACCTCGGTTGAAACTGACAGGTGGTTTGTTACGCATGCTAATGCAAAGGAGCCTATATACCTTTGGCTCGGCTGCTGTAACAGGGAATATAAAGGGCAGCATAATTTAGGAGTTTAGTGAACTTGCAACATTTACTATTTTCCCTTCTTACGTAAATATTTTTCTTTTTAATTCTAAATCAATCTTTTTCAATTTTTTGTTTGTATTCTTTTCTTGCTTAAATCTATAACTACAAAAAACACATACATAAACTAAAA'
    s1 = sc.intergenic_sequence("YDR051C", "TPI1")
    s2 = sc.promoter("TPI1")
    assert tp.lower() == str(s1.seq).lower() == str(s2.seq).lower()

def test_TPI1_genbank_codes():
    assert sc.intergenic_sequence_genbank_accession("YDR051C", "TPI1") == 'BK006938.2 REGION: complement(556473..557055)'
    assert sc.promoter_genbank_accession("TPI1") == 'BK006938.2 REGION: complement(556473..557055)'

def test_fun26():
    assert len(str(sc.locus("fun26").seq)) == 3554
    assert len(str(sc.cds("fun26").seq))   == 1554
    assert sc.upstream_gene("fun26") == 'YAL021C'
    assert sc.systematic_name("fun26") == 'YAL022C'
    assert sc.downstream_gene("fun26") == 'YAL023C'
    assert str(sc.intergenic_sequence("YAL021C","YAL023C").seq).lower() in str(sc.locus("fun26").seq).lower()
    assert str(  sc.promoter("fun26").seq) in str(sc.locus("fun26").seq)
    assert str(  sc.terminator("fun26").seq) in str(sc.locus("fun26").seq)
    assert str(  sc.promoter("pmt2").seq) in str(sc.locus("pmt2").seq)
    assert str(  sc.terminator("lte1").seq) in str(sc.locus("lte1").seq)


def test_promoter_promoter():
    assert str(sc.promoter("dep1").seq) == str(sc.promoter("syn8").seq.reverse_complement())
    assert str(sc.promoter("spo7").seq) == str(sc.promoter("mdm10").seq.reverse_complement())

def test_teminator_terminator():
    assert str(sc.terminator("fun14").seq) == str(sc.terminator("erp2").seq.reverse_complement())

def test_promoter_terminator():
    assert str(sc.promoter("cys3").seq) == str(sc.terminator("dep1").seq)

def test_terminator_promoter():
    assert str(sc.promoter("ccr4").seq) == str(sc.terminator("ats1").seq)

def test_tandem_bidirectional():
    assert sc.tandem("tpi1") == (not sc.bidirectional("tpi1"))
    assert sc.tandem("gal1") == (not sc.bidirectional("gal1"))

def test_misc():
    assert str(sc.promoter("CLN3").seq) in str(sc.locus("CLN3",2000,2000).seq)
    assert str(sc.terminator("CLN3").seq) in str(sc.locus("CLN3",2000,2000).seq)
    assert str(sc.terminator("cyc3").seq) in str(sc.locus("cyc3",2000,2000).seq)
    assert str(sc.promoter("cyc3").seq) in str(sc.locus("cyc3",2500,2500).seq)
    assert str(sc.terminator("cyc3").seq) == str(sc.promoter("CLN3").seq)
    assert str(sc.promoter("jen1").seq) == str(sc.promoter("sry1").seq.reverse_complement())
    assert str(sc.promoter("osm1").seq) == str(sc.terminator("isy1").seq)
    assert str(sc.terminator("cyc1").seq) in str(sc.locus("cyc1").seq)
    assert str(sc.terminator("utr1").seq) in str(sc.locus("utr1").seq)
    assert str(  sc.promoter("cdc24").seq) in str(sc.locus("cdc24").seq)
    assert str(sc.terminator("cdc24").seq) in str(sc.locus("cdc24").seq)

if __name__ == '__main__':
    nose.runmodule()
