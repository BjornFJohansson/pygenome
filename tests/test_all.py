#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
test pygenome
'''
import pytest
from pygenome import sg

def test_TEF1():
    tp = 'ACAATGCATACTTTGTACGTTCAAAATACAATGCAGTAGATATATTTATGCATATTACATATAATACATATCACATAGGAAGCAACAGGCGCGTTGGACTTTTAATTTTCGAGGACCGCGAATCCTTACATCACACCCAATCCCCCACAAGTGATCCCCCACACACCATAGCTTCAAAATGTTTCTACTCCTTTTTTACTCTTCCAGATTTTCTCGGACTCCGCGCATCGCCGTACCACTTCAAAACACCCAAGCACAGCATACTAAATTTCCCCTCTTTCTTCCTCTAGGGTGTCGTTAATTACCCGTACTAAAGGTTTGGAAAAGAAAAAAGAGACCGCCTCGTTTCTTTTTCTTCGTCGAAAAAGGCAATAAAAATTTTTATCACGTTTCTTTTTCTTGAAAATTTTTTTTTTTGATTTTTTTCTCTTTCGATGACCTCCCATTGATATTTAAGTTAATAAACGGTCTTCAATTTCTCAAGTTTCAGTTTCATTTTTCTTGTTCTATTACAACTTTTTTTACTTCTTGCTCATTAGAAAGAAAGCATAGCAATCTAATCTAAGTTTTAATTACAAA'
    s1 = sg.gene["YPR079W"].terminator
    s2 = sg.gene["TEF1"].promoter
    s3 = sg.gene["YPR080W"].promoter
    assert tp.lower() == str(s1.seq).lower() == str(s2.seq).lower() == str(s3.seq).lower()

def test_TPI1():
    tp = 'TGTTTAAAGATTACGGATATTTAACTTACTTAGAATAATGCCATTTTTTTGAGTTATAATAATCCTACGTTAGTGTGAGCGGGATTTAAACTGTGAGGACCTTAATACATTCAGACACTTCTGCGGTATCACCCTACTTATTCCCTTCGAGATTATATCTAGGAACCCATCAGGTTGGTGGAAGATTACCCGTTCTAAGACTTTTCAGCTTCCTCTATTGATGTTACACCTGGACACCCCTTTTCTGGCATCCAGTTTTTAATCTTCAGTGGCATGTGAGATTCTCCGAAATTAATTAAAGCAATCACACAATTCTCTCGGATACCACCTCGGTTGAAACTGACAGGTGGTTTGTTACGCATGCTAATGCAAAGGAGCCTATATACCTTTGGCTCGGCTGCTGTAACAGGGAATATAAAGGGCAGCATAATTTAGGAGTTTAGTGAACTTGCAACATTTACTATTTTCCCTTCTTACGTAAATATTTTTCTTTTTAATTCTAAATCAATCTTTTTCAATTTTTTGTTTGTATTCTTTTCTTGCTTAAATCTATAACTACAAAAAACACATACATAAACTAAAA'
    s1 = sg.gene["YDR051C"].terminator
    s2 = sg.gene["TPI1"].promoter
    s3 = sg.gene["YDR050C"].promoter
    assert tp.lower() == str(s1.seq).lower() == str(s2.seq).lower() == str(s3.seq).lower()

def test_TEF1_genbank_accession():
    assert sg.gene["TEF1"].promoter.description == 'BK006949.2 REGION: 700015..700593'

def test_TPI1_genbank_accession():
    assert sg.gene["TPI1"].promoter.description == 'BK006938.2 REGION: complement(556473..557055)'

def test_fun26():
    assert len(sg.gene["FUN26"].locus()) == 3554
    assert len(sg.gene["FUN26"].cds)   == 1554
    assert sg.gene["FUN26"].sys == 'YAL022C'
    assert sg.gene["FUN26"].upstream_gene.sys == 'YAL021C'
    assert sg.gene["FUN26"].upstream_gene.sys == 'YAL021C'
    assert sg.gene["FUN26"].downstream_gene.sys == 'YAL023C'

    assert str(  sg.gene["FUN26"].promoter.seq)   in str(sg.gene["FUN26"].locus().seq)
    assert str(  sg.gene["FUN26"].terminator.seq) in str(sg.gene["FUN26"].locus().seq)
    assert str(  sg.gene["PMT2"].promoter.seq)    in str(sg.gene["PMT2"].locus().seq)
    assert str(  sg.gene["LTE1"].terminator.seq)  in str(sg.gene["LTE1"].locus().seq)



def test_promoter_promoter():
    assert str(sg.gene["DEP1"].promoter.seq) == str(sg.gene["SYN8"].promoter.seq.reverse_complement())
    assert str(sg.gene["SPO7"].promoter.seq) == str(sg.gene["MDM10"].promoter.seq.reverse_complement())

def test_teminator_terminator():
    assert str(sg.gene["FUN14"].terminator.seq) == str(sg.gene["ERP2"].terminator.seq.reverse_complement())

def test_promoter_terminator():
    assert str(sg.gene["CYS3"].promoter.seq) == str(sg.gene["DEP1"].terminator.seq)

def test_terminator_promoter():
    assert str(sg.gene["CCR4"].promoter.seq) == str(sg.gene["ATS1"].terminator.seq)

def test_tandem_bidirectional():
    assert sg.gene["TPI1"].tandem == (not sg.gene["TPI1"].bidirectional)
    assert sg.gene["GAL1"].tandem == (not sg.gene["GAL1"].bidirectional)

def test_misc():

    assert str(sg.gene["CLN3"].promoter.seq) in str(sg.gene["CLN3"].locus(2000, 2000).seq)
    assert str(sg.gene["CLN3"].terminator.seq) in str(sg.gene["CLN3"].locus(2000, 2000).seq)

    assert str(sg.gene["CYC3"].terminator.seq) in str(sg.gene["CYC3"].locus(2000, 2000).seq)

    assert str(sg.gene["CYC3"].promoter.seq)    in str(sg.gene["CYC3"].locus(2500,2500).seq)
    assert str(sg.gene["CYC3"].terminator.seq)  == str(sg.gene["CLN3"].promoter.seq)
    assert str(sg.gene["JEN1"].promoter.seq)    == str(sg.gene["SRY1"].promoter.seq.reverse_complement())
    assert str(sg.gene["OSM1"].promoter.seq)    == str(sg.gene["ISY1"].terminator.seq)
    assert str(sg.gene["CYC1"].terminator.seq)  in str(sg.gene["CYC1"].locus().seq)
    assert str(sg.gene["UTR1"].terminator.seq)  in str(sg.gene["UTR1"].locus().seq)
    assert str(sg.gene["CDC24"].promoter.seq)   in str(sg.gene["CDC24"].locus().seq)
    assert str(sg.gene["CDC24"].terminator.seq) in str(sg.gene["CDC24"].locus().seq)

    assert str( sg.intergenic_sequence("YAL021C","YAL023C").seq).lower() in str(sg.gene["FUN26"].locus().seq).lower()

def test_kanmx4():

    s = sg.gene["CYC1"].deletion_locus



    text = '''
    gaggcaccagcgtcagcattttcaaaggtgtgttcttcgtcagacatgttttagtgtgtgaatgaaataggtgtatgttttctttttgctagacaataattaggaacaaggtaagggaactaaagtgtagaataagattaaaaaagaagaacaagttgaaaaggcaagttgaaatttcaagaaaaaagtcaattgaagtacagtaaattgacctgaatatatctgagttccgacaacaatgagtttaccaaagagaacaatggaataggaaactttgaacgaagaaaggaaagcaggaaaggaaaaaatttttaggctcgagaacaatagggcgaaaaaacaggcaacgaacgaacaatggaaaaacgaaaaaaaaaaaaaaaaacacagaaaagaatgcagaaagatgtcaactgaaaaaaaaaaaggtgaacacaggaaaaaaaataaaaaaaaaaaaaaaaaaaggaggacgaaacaaaaaagtgaaaaaaaatgaaaatttttttggaaaaccaagaaatgaattatatttccgtgtgagacgacatcgtcgaatatgattcagggtaacagtattgatgtaatcaatttcctacctgaatctaaaattcccgggagcaagatcaagatgttttcaccgatctttccggtctctttggccggggtttacggacgatggcagaagaccaaagcgccagttcatttggcgagcgttggttggtggatcaagcccacgcgtaggcaatcctcgagcagatccgccaggcgtgtatatatagcgtggatggccaggcaactttagtgctgacacatacaggcatatatatatgtgtgcgacgacacatgatcatatggcatgcatgtgctctgtatgtatataaaactcttgttttcttcttttctctaaatattctttccttatacattaggacctttgcagcataaattactatac



    TTCTATAGACACACAAACACAAATACA
                               CACACTAAATTAATAATGGATGTCCACGAGGTCTCTATATCGGGATCAGCCTGCCTCGTACGCTGCAGGTCGAC
                                                                                       GGAT
    CCCCGGGTTAATTAAGGCGCGCCAGATCTGTTTAGCTTGCCTCGTCCCCGCCGGGTCACCCGGCCAGCGACATGGAGGCC
    CAGAATACCCTCCTTGACAGTCTTGACGTGCGCAGCTCAGGGGCATGATGTGACTGTCGCCCGTACATTTAGCCCATACATCCCCATGTATAATCATTTGCATCC
    ATACATTTTGATGGCCGCACGGCGCGAAGCAAAAATTACGGCTCCTCGCTGCAGACCTGCGAGCAGGGAAACGCTCCCCTCACAGACGCGTTGAATTGTCCCCAC
    GCCGCGCCCCTGTAGAGAAATATAAAAGGTTAGGATTTGCCACTGAGGTTCTTCTTTCATATACTTCCTTTTAAAATCTTGCTAGGATACAGTTCTCACATCACAT
    CCGAACATAAACAACCATGGGTAAGGAAAAGACTCACGTTTCGAGGCCGCGATTAAATTCCAACATGGATGCTGATTTATATGGGTATAAATGGGCTCGCGATAATG
    TCGGGCAATCAGGTGCGACAATCTATCGATTGTATGGGAAGCCCGATGCGCCAGAGTTGTTTCTGAAACATGGCAAAGGTAGCGTTGCCAATGATGTTACAGATG
    AGATGGTCAGACTAAACTGGCTGACGGAATTTATGCCTCTTCCGACCATCAAGCATTTTATCCGTACTCCTGATGATGCATGGTTACTCACCACTGCGATCCCCG
    GCAAAACAGCATTCCAGGTATTAGAAGAATATCCTGATTCAGGTGAAAATATTGTTGATGCGCTGGCAGTGTTCCTGCGCCGGTTGCATTCGATTCCTGTTTGTAA
    TTGTCCTTTTAACAGCGATCGCGTATTTCGTCTCGCTCAGGCGCAATCACGAATGAATAACGGTTTGGTTGATGCGAGTGATTTTGATGACGAGCGTAATGGCTGGC
    CTGTTGAACAAGTCTGGAAAGAAATGCATAAGCTTTTGCCATTCTCACCGGATTCAGTCGTCACTCATGGTGATTTCTCACTTGATAACCTTATTTTTGACGAGGGG
    AAATTAATAGGTTGTATTGATGTTGGACGAGTCGGAATCGCAGACCGATACCAGGATCTTGCCATCCTATGGAACTGCCTCGGTGAGTTTTCTCCTTCATTACAGAA
    ACGGCTTTTTCAAAAATATGGTATTGATAATCCTGATATGAATAAATTGCAGTTTCATTTGATGCTCGATGAGTTTTTCTAATCAGTACTGACAATAAAAAGATTCT
    TGTTTTCAAGAACTTGTCATTTGTATAGTTTTTTTATATTGTAGTTGTTCTATTTTAATCAAATGTTAGCGTGATTTATATTTTTTTTCGCCTCGACATCATCTGCCC
    AGATGCGAAGTTAAGTGCGCAGAAAGTAATATCATGCGTCAATCGTATGTGAATGCTGGTCGCTATACTGCTGTCGATTCGATACTAACGCCGCCATCCAGTGT
    CGAAAA
          CGAGCTCGAATTCATCGATGAGCATCCTTATAGCCTCTTCTACGAGACCGACACCG
                                                                  TAAACAGGCCCCTTTTCCTTTGTCGATATCATGTAATTAGTTATG

    tcacgcttacattcacgccctcctcccacatccgctctaaccgaaaaggaaggagttagacaacctgaagtctaggtccctatttattttttttaatagttatgttagtattaagaacgttatttatatttcaaatttttcttttttttctgtacaaacgcgtgtacgcatgtaacattatactgaaaaccttgcttgagaaggttttgggacgctcgaaggctttaatttgcaagcttcgcagtttacactctcatcgtcgctctcatcatcgcttccgttgttgttttccttagtagcgtctgcttccagagagtatttatctcttattacctctaaaggttctgcttgatttctgactttgttcgcctcatgtgcatatttttcttggttcttttgggacaaaatatgcgtaaaggacttttgttgttccctcacattccagtttagttgtcgactgatactgttaataaactcatcgggcgaggcttccacggttggaaaagcatatgggctggcgcatatggttataaaatcacctttttgcaattcaattctatctttcccatcaaaagccgcccatgctggagcccttgacttcatcgagactttcacttttaaatttatactttctggtaagatgatgggtctgaaactcaatgcatgtggacaaatgggtgttaaagcgattgcattgacggttgggcataccaatgacccacctgcactcaaagaataggccgtggacccagtcggagtagcagcaatcagtccgtccgcctgcgcaacggtcattaatgagccgtcaccatacaattctaacatggatagaaaaggacttggaccacgatcgatggtcacttcgttcaaaatgtggtgtgtgcttagtttttccaccacacatattttcttccccgtgtttgggtctacttcagggcggtgtctacgataaattgtg

    '''

    text = "".join([c.strip() for c in text])

    assert text.lower() == str(s.seq).lower()


if __name__ == '__main__':
    pytest.cmdline.main([__file__, "-v", "-s"])


# CYC1 deletion locus with the four primers used to make the cassette
'''

gaggcaccagcgtcagcattttcaaaggtgtgttcttcgtcagacatgttttagtgtgtgaatgaaataggtgtatgttttctttttgctagacaataat
taggaacaaggtaagggaactaaagtgtagaataagattaaaaaagaagaacaagttgaaaaggcaagttgaaatttcaagaaaaaagtcaattgaagta
cagtaaattgacctgaatatatctgagttccgacaacaatgagtttaccaaagagaacaatggaataggaaactttgaacgaagaaaggaaagcaggaaa
ggaaaaaatttttaggctcgagaacaatagggcgaaaaaacaggcaacgaacgaacaatggaaaaacgaaaaaaaaaaaaaaaaacacagaaaagaatgca
gaaagatgtcaactgaaaaaaaaaaaggtgaacacaggaaaaaaaataaaaaaaaaaaaaaaaaaaggaggacgaaacaaaaaagtgaaaaaaaatgaaaat
ttttttggaaaaccaagaaatgaattatatttccgtgtgagacgacatcgtcgaatatgattcagggtaacagtattgatgtaatcaatttcctacctgaa
tctaaaattcccgggagcaagatcaagatgttttcaccgatctttccggtctctttggccggggtttacggacgatggcagaagaccaaagcgccagttcat
ttggcgagcgttggttggtggatcaagcccacgcgtaggcaatcctcgagcagatccgccaggcgtgtatatatagcgtggatggccaggcaactttagtg
ctgacacatacaggcatatatatatgtgtgcgacgacacatgatcatatggcatgcatgtgctctgtatgtatataaaactcttgttttcttcttttctct
aaatattctttccttatacattaggacctttgcagcataaattactatac



TTCTATAGACACACAAACACAAATACACACACTAAATTAATAATG
                           CACACTAAATTAATAATGGATGTCCACGAGGTCTCTATATCGGGATCAGCCTGCCTCGTACGCTGCAGGTCGAC
                                                                                   CGTACGCTGCAGGTCGACGGAT
CCCCGGGTTAATTAAGGCGCGCCAGATCTGTTTAGCTTGCCTCGTCCCCGCCGGGTCACCCGGCCAGCGACATGGAGGCC
CAGAATACCCTCCTTGACAGTCTTGACGTGCGCAGCTCAGGGGCATGATGTGACTGTCGCCCGTACATTTAGCCCATACATCCCCATGTATAATCATTTGCATCC
ATACATTTTGATGGCCGCACGGCGCGAAGCAAAAATTACGGCTCCTCGCTGCAGACCTGCGAGCAGGGAAACGCTCCCCTCACAGACGCGTTGAATTGTCCCCAC
GCCGCGCCCCTGTAGAGAAATATAAAAGGTTAGGATTTGCCACTGAGGTTCTTCTTTCATATACTTCCTTTTAAAATCTTGCTAGGATACAGTTCTCACATCACAT
CCGAACATAAACAACCATGGGTAAGGAAAAGACTCACGTTTCGAGGCCGCGATTAAATTCCAACATGGATGCTGATTTATATGGGTATAAATGGGCTCGCGATAATG
TCGGGCAATCAGGTGCGACAATCTATCGATTGTATGGGAAGCCCGATGCGCCAGAGTTGTTTCTGAAACATGGCAAAGGTAGCGTTGCCAATGATGTTACAGATG
AGATGGTCAGACTAAACTGGCTGACGGAATTTATGCCTCTTCCGACCATCAAGCATTTTATCCGTACTCCTGATGATGCATGGTTACTCACCACTGCGATCCCCG
GCAAAACAGCATTCCAGGTATTAGAAGAATATCCTGATTCAGGTGAAAATATTGTTGATGCGCTGGCAGTGTTCCTGCGCCGGTTGCATTCGATTCCTGTTTGTAA
TTGTCCTTTTAACAGCGATCGCGTATTTCGTCTCGCTCAGGCGCAATCACGAATGAATAACGGTTTGGTTGATGCGAGTGATTTTGATGACGAGCGTAATGGCTGGC
CTGTTGAACAAGTCTGGAAAGAAATGCATAAGCTTTTGCCATTCTCACCGGATTCAGTCGTCACTCATGGTGATTTCTCACTTGATAACCTTATTTTTGACGAGGGG
AAATTAATAGGTTGTATTGATGTTGGACGAGTCGGAATCGCAGACCGATACCAGGATCTTGCCATCCTATGGAACTGCCTCGGTGAGTTTTCTCCTTCATTACAGAA
ACGGCTTTTTCAAAAATATGGTATTGATAATCCTGATATGAATAAATTGCAGTTTCATTTGATGCTCGATGAGTTTTTCTAATCAGTACTGACAATAAAAAGATTCT
TGTTTTCAAGAACTTGTCATTTGTATAGTTTTTTTATATTGTAGTTGTTCTATTTTAATCAAATGTTAGCGTGATTTATATTTTTTTTCGCCTCGACATCATCTGCCC
AGATGCGAAGTTAAGTGCGCAGAAAGTAATATCATGCGTCAATCGTATGTGAATGCTGGTCGCTATACTGCTGTCGATTCGATACTAACGCCGCCATCCAGTGT
CGAAAACGAGCTCGAATTCATCGAT
      CGAGCTCGAATTCATCGATGAGCATCCTTATAGCCTCTTCTACGAGACCGACACCGTAAACAGGCCCCTTTTCC
                                                              TAAACAGGCCCCTTTTCCTTTGTCGATATCATGTAATTAGTTATG

tcacgcttacattcacgccctcctcccacatccgctctaaccgaaaaggaaggagttagacaacctgaagtctaggtccctatttattttttttaatagttatgttagt
attaagaacgttatttatatttcaaatttttcttttttttctgtacaaacgcgtgtacgcatgtaacattatactgaaaaccttgcttgagaaggttttgggacgctcg
aaggctttaatttgcaagcttcgcagtttacactctcatcgtcgctctcatcatcgcttccgttgttgttttccttagtagcgtctgcttccagagagtatttatctct
tattacctctaaaggttctgcttgatttctgactttgttcgcctcatgtgcatatttttcttggttcttttgggacaaaatatgcgtaaaggacttttgttgttccct
cacattccagtttagttgtcgactgatactgttaataaactcatcgggcgaggcttccacggttggaaaagcatatgggctggcgcatatggttataaaatcaccttt
ttgcaattcaattctatctttcccatcaaaagccgcccatgctggagcccttgacttcatcgagactttcacttttaaatttatactttctggtaagatgatgggtct
gaaactcaatgcatgtggacaaatgggtgttaaagcgattgcattgacggttgggcataccaatgacccacctgcactcaaagaataggccgtggacccagtcggagt
agcagcaatcagtccgtccgcctgcgcaacggtcattaatgagccgtcaccatacaattctaacatggatagaaaaggacttggaccacgatcgatggtcacttcgtt
caaaatgtggtgtgtgcttagtttttccaccacacatattttcttccccgtgtttgggtctacttcagggcggtgtctacgataaattgtg

'''