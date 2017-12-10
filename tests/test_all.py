#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
test pygenome
'''
import pytest
import os
from pygenome import sg
from pygenome._pretty import pretty_str
from pygenome.systematic_name import _systematic_name


from pygenome._data import  _data_files, _data_urls
from pygenome.update import updater
import requests_mock as rm_module
import shutil, io, pathlib

@pytest.fixture
def requests_mock(request):
    m = rm_module.Mocker()
    m.start()
    request.addfinalizer(m.stop)
    return m



def test_update(requests_mock):
    data_dir = os.path.join(os.getenv("pygenome_data_dir"), "Saccharomyces_cerevisiae")
    tmp_data_dir = os.path.join(os.getenv("pygenome_data_dir"), "temp") # copy all data files to another directory
    shutil.copytree(data_dir, tmp_data_dir)
    
    for fn, url in zip(_data_files,_data_urls): 
        # set data file to be really old....
        path = pathlib.Path( os.path.join(data_dir, fn) )
        os.utime(str(path), times=(path.stat().st_atime, 946684800))  #946684800 = Saturday 1st January 2000 12:00:00 AM
        flo = io.BytesIO(b"some text data")  # These files will be deemed newer and downloaded
        requests_mock.get(url, 
                          headers={'last-modified'  : 'Mon, 01 Jan 2001 00:00:00 GMT', #978307200
                                   'content-length' : "100"}, 
                          body = flo)
    updater()  # should "download" the newer files
    for fn, url in zip(_data_files,_data_urls): 
        # set data file to be the same age remote
        path = pathlib.Path( os.path.join(data_dir, fn) )
        os.utime(str(path), times=(path.stat().st_atime, 978307200))
        flo = io.BytesIO(b"some text data that will not be used")
        requests_mock.get(url, 
                          headers={'last-modified'  : 'Mon, 01 Jan 2001 00:00:00 GMT', #978307200
                                   'content-length' : "100"}, 
                          body = flo)
    updater() # should keep the local files  
    for fn, url in zip(_data_files,_data_urls):
        path = pathlib.Path( os.path.join(data_dir, fn) )
        path.unlink()
        flo = io.BytesIO(b"some text data that will be used")
        requests_mock.get(url, 
                          headers={'last-modified'  : 'Mon, 01 Jan 2001 00:00:00 GMT', #978307200
                                   'content-length' : "100"}, 
                          body = flo)
        
    updater()    
    shutil.rmtree(data_dir)
    shutil.copytree(tmp_data_dir, data_dir)

def test_pretty():
    from unittest.mock import MagicMock
    pp = MagicMock()
    x=pretty_str("abc")
    x._repr_pretty_(pp, None)
    pp.text.assert_any_call("abc")
    
def test_names():
    x=_systematic_name("TDH3")
    with pytest.raises(KeyError):
        _systematic_name("NOGENE")


def test_TEF1():
    tp = 'ACAATGCATACTTTGTACGTTCAAAATACAATGCAGTAGATATATTTATGCATATTACATATAATACATATCACATAGGAAGCAACAGGCGCGTTGGACTTTTAATTTTCGAGGACCGCGAATCCTTACATCACACCCAATCCCCCACAAGTGATCCCCCACACACCATAGCTTCAAAATGTTTCTACTCCTTTTTTACTCTTCCAGATTTTCTCGGACTCCGCGCATCGCCGTACCACTTCAAAACACCCAAGCACAGCATACTAAATTTCCCCTCTTTCTTCCTCTAGGGTGTCGTTAATTACCCGTACTAAAGGTTTGGAAAAGAAAAAAGAGACCGCCTCGTTTCTTTTTCTTCGTCGAAAAAGGCAATAAAAATTTTTATCACGTTTCTTTTTCTTGAAAATTTTTTTTTTTGATTTTTTTCTCTTTCGATGACCTCCCATTGATATTTAAGTTAATAAACGGTCTTCAATTTCTCAAGTTTCAGTTTCATTTTTCTTGTTCTATTACAACTTTTTTTACTTCTTGCTCATTAGAAAGAAAGCATAGCAATCTAATCTAAGTTTTAATTACAAA'
    s1 = sg.sysgene["YPR079W"].terminator
    s2 = sg.stdgene["TEF1"].promoter
    s3 = sg.sysgene["YPR080W"].promoter
    assert tp.lower() == str(s1.seq).lower() == str(s2.seq).lower() == str(s3.seq).lower()

def test_TPI1():
    tp = 'TGTTTAAAGATTACGGATATTTAACTTACTTAGAATAATGCCATTTTTTTGAGTTATAATAATCCTACGTTAGTGTGAGCGGGATTTAAACTGTGAGGACCTTAATACATTCAGACACTTCTGCGGTATCACCCTACTTATTCCCTTCGAGATTATATCTAGGAACCCATCAGGTTGGTGGAAGATTACCCGTTCTAAGACTTTTCAGCTTCCTCTATTGATGTTACACCTGGACACCCCTTTTCTGGCATCCAGTTTTTAATCTTCAGTGGCATGTGAGATTCTCCGAAATTAATTAAAGCAATCACACAATTCTCTCGGATACCACCTCGGTTGAAACTGACAGGTGGTTTGTTACGCATGCTAATGCAAAGGAGCCTATATACCTTTGGCTCGGCTGCTGTAACAGGGAATATAAAGGGCAGCATAATTTAGGAGTTTAGTGAACTTGCAACATTTACTATTTTCCCTTCTTACGTAAATATTTTTCTTTTTAATTCTAAATCAATCTTTTTCAATTTTTTGTTTGTATTCTTTTCTTGCTTAAATCTATAACTACAAAAAACACATACATAAACTAAAA'
    s1 = sg.sysgene["YDR051C"].terminator
    s2 = sg.stdgene["TPI1"].promoter
    s3 = sg.sysgene["YDR050C"].promoter
    assert tp.lower() == str(s1.seq).lower() == str(s2.seq).lower() == str(s3.seq).lower()

def test_TEF1_genbank_accession():
    assert sg.stdgene["TEF1"].promoter.description == 'BK006949.2 REGION: 700015..700593'

def test_TPI1_genbank_accession():
    assert sg.stdgene["TPI1"].promoter.description == 'BK006938.2 REGION: complement(556473..557055)'

def test_fun26():
    assert len(sg.stdgene["FUN26"].locus()) == 3554
    assert len(sg.stdgene["FUN26"].cds)   == 1554
    assert sg.stdgene["FUN26"].sys == 'YAL022C'
    assert sg.stdgene["FUN26"].upstream_gene.sys == 'YAL021C'
    assert sg.stdgene["FUN26"].upstream_gene.sys == 'YAL021C'
    assert sg.stdgene["FUN26"].downstream_gene.sys == 'YAL023C'

    assert str(  sg.stdgene["FUN26"].promoter.seq)   in str(sg.stdgene["FUN26"].locus().seq)
    assert str(  sg.stdgene["FUN26"].terminator.seq) in str(sg.stdgene["FUN26"].locus().seq)
    assert str(  sg.stdgene["PMT2"].promoter.seq)    in str(sg.stdgene["PMT2"].locus().seq)
    assert str(  sg.stdgene["LTE1"].terminator.seq)  in str(sg.stdgene["LTE1"].locus().seq)



def test_promoter_promoter():
    assert str(sg.stdgene["DEP1"].promoter.seq) == str(sg.stdgene["SYN8"].promoter.seq.reverse_complement())
    assert str(sg.stdgene["SPO7"].promoter.seq) == str(sg.stdgene["MDM10"].promoter.seq.reverse_complement())

def test_terminator_terminator():
    assert str(sg.stdgene["FUN14"].terminator.seq) == str(sg.stdgene["ERP2"].terminator.seq.reverse_complement())

def test_promoter_terminator():
    assert str(sg.stdgene["CYS3"].promoter.seq) == str(sg.stdgene["DEP1"].terminator.seq)

def test_terminator_promoter():
    assert str(sg.stdgene["CCR4"].promoter.seq) == str(sg.stdgene["ATS1"].terminator.seq)

def test_tandem_bidirectional():
    assert sg.stdgene["TPI1"].tandem == (not sg.stdgene["TPI1"].bidirectional)
    assert sg.stdgene["GAL1"].tandem == (not sg.stdgene["GAL1"].bidirectional)

def test_misc():

    assert str(sg.stdgene["CLN3"].promoter.seq) in str(sg.stdgene["CLN3"].locus(2000, 2000).seq)
    assert str(sg.stdgene["CLN3"].terminator.seq) in str(sg.stdgene["CLN3"].locus(2000, 2000).seq)

    assert str(sg.stdgene["CYC3"].terminator.seq) in str(sg.stdgene["CYC3"].locus(2000, 2000).seq)

    assert str(sg.stdgene["CYC3"].promoter.seq)    in str(sg.stdgene["CYC3"].locus(2500,2500).seq)
    assert str(sg.stdgene["CYC3"].terminator.seq)  == str(sg.stdgene["CLN3"].promoter.seq)
    assert str(sg.stdgene["JEN1"].promoter.seq)    == str(sg.stdgene["SRY1"].promoter.seq.reverse_complement())
    assert str(sg.stdgene["OSM1"].promoter.seq)    == str(sg.stdgene["ISY1"].terminator.seq)
    assert str(sg.stdgene["CYC1"].terminator.seq)  in str(sg.stdgene["CYC1"].locus().seq)
    assert str(sg.stdgene["UTR1"].terminator.seq)  in str(sg.stdgene["UTR1"].locus().seq)
    assert str(sg.stdgene["CDC24"].promoter.seq)   in str(sg.stdgene["CDC24"].locus().seq)
    assert str(sg.stdgene["CDC24"].terminator.seq) in str(sg.stdgene["CDC24"].locus().seq)

    from pygenome.intergenic import intergenic_sequence
    intseq = intergenic_sequence("YAL021C","YAL023C")
    assert str( intseq.seq).lower() in str(sg.stdgene["FUN26"].locus().seq).lower()
    

    with pytest.raises(Exception):
        intergenic_sequence("YAL021C","YBL023C")
    

def test_kanmx4():

    s = sg.stdgene["CYC1"].deletion_locus



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
    
def test_pickle():
    data_dir = os.path.join(os.getenv("pygenome_data_dir"), "Saccharomyces_cerevisiae")
    import pathlib
    for p in pathlib.Path(data_dir).glob("*.pickle"): p.unlink()
    
    from pygenome._pickle_primers import pickle_primers
    pickle_primers()

    from pygenome._pickle_lists   import _pickle_lists
    _pickle_lists()
    
    from pygenome._pickle_primers import pickle_orfs_not_deleted
    pickle_orfs_not_deleted()

    from pygenome._pickle_genes   import _pickle_genes
    _pickle_genes()
    

