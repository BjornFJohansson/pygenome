#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
test pygenome
'''
import pytest
from pathlib import Path
import requests_mock as rm_module
from tempfile import TemporaryDirectory
import shutil
import io
import os
from urllib.parse import urlparse as up
from os.path import basename as bn
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

#  https://realpython.com/python-mock-library


@pytest.fixture
def requests_mock(request):
    m = rm_module.Mocker()
    m.start()
    request.addfinalizer(m.stop)
    return m


def test_check_and_download(requests_mock, monkeypatch):

    def Bio_SeqIO_read(pth, fmat):
        result = SeqRecord(Seq("GATC"))
        return result   

    monkeypatch.setattr("Bio.SeqIO.read",
                        Bio_SeqIO_read,
                        raising=True)
    
    from pygenome.saccharomyces_cerevisiae.S288C import check_data_files
    from pygenome.saccharomyces_cerevisiae.S288C import update_data_files
    from pygenome.saccharomyces_cerevisiae.S288C import data_dir
    from pygenome.saccharomyces_cerevisiae.S288C import chromosome_urls

    tmp_data_dir = Path(TemporaryDirectory().name) # tempdir for data files
    
    shutil.copytree(data_dir, tmp_data_dir)
    
    # set local data file to be really old.... 
    # Saturday 1st January 2000 12:00:00 AM   =  946684800
    # new files should be downloaded

    chromosome_urls = chromosome_urls.splitlines()

    for url in chromosome_urls:

        fn = bn(up(url).path)
        gb_pth = data_dir/fn
        fa_pth = gb_pth.with_suffix(".fasta")
        fe_pth = fa_pth.with_name(f"{fa_pth.stem}_feature_tuple.pickle")

        os.utime(fa_pth, times=(fa_pth.stat().st_atime, 946684800))

        flo = io.BytesIO(b"some text data") # mock remote file is newer
        requests_mock.get(url,
                          headers={  # 978307200
                           'last-modified': 'Mon, 01 Jan 2001 00:00:00 GMT',
                           'content-length': "100"},
                          body=flo)

    check_data_files()
    
    update_data_files()

    shutil.rmtree(data_dir)
    shutil.copytree(tmp_data_dir, data_dir)
    shutil.rmtree(tmp_data_dir)


def test_Gene():
    from pygenome.saccharomyces_cerevisiae.S288C import gene_dicts
    genes, stdgenes = gene_dicts()

    g = genes["YGR192C"]
    s = stdgenes["TDH3"]
    assert g == s

    with pytest.raises(KeyError):
        assert genes["NOGENE"]
    with pytest.raises(KeyError):
        assert stdgenes["NOGENE"]

    tp = ("ACAATGCATACTTTGTACGTTCAAAATACAATGCAGTAGATATATTTATGCATATTACATA"
          "TAATACATATCACATAGGAAGCAACAGGCGCGTTGGACTTTTAATTTTCGAGGACCGCGAA"
          "TCCTTACATCACACCCAATCCCCCACAAGTGATCCCCCACACACCATAGCTTCAAAATGTT"
          "TCTACTCCTTTTTTACTCTTCCAGATTTTCTCGGACTCCGCGCATCGCCGTACCACTTCAA"
          "AACACCCAAGCACAGCATACTAAATTTCCCCTCTTTCTTCCTCTAGGGTGTCGTTAATTAC"
          "CCGTACTAAAGGTTTGGAAAAGAAAAAAGAGACCGCCTCGTTTCTTTTTCTTCGTCGAAAA"
          "AGGCAATAAAAATTTTTATCACGTTTCTTTTTCTTGAAAATTTTTTTTTTTGATTTTTTTC"
          "TCTTTCGATGACCTCCCATTGATATTTAAGTTAATAAACGGTCTTCAATTTCTCAAGTTTC"
          "AGTTTCATTTTTCTTGTTCTATTACAACTTTTTTTACTTCTTGCTCATTAGAAAGAAAGCA"
          "TAGCAATCTAATCTAAGTTTTAATTACAAA")
    s1 = genes["YPR079W"].terminator
    s2 = stdgenes["TEF1"].promoter
    s3 = genes["YPR080W"].promoter
    assert tp.lower() == str(s1.seq).lower()
    assert str(s1.seq).lower() == str(s2.seq).lower()
    assert str(s2.seq).lower() == str(s3.seq).lower()

    tp = ("TGTTTAAAGATTACGGATATTTAACTTACTTAGAATAATGCCATTTTTTTGAGTTATAATA"
          "ATCCTACGTTAGTGTGAGCGGGATTTAAACTGTGAGGACCTTAATACATTCAGACACTTCT"
          "GCGGTATCACCCTACTTATTCCCTTCGAGATTATATCTAGGAACCCATCAGGTTGGTGGAA"
          "GATTACCCGTTCTAAGACTTTTCAGCTTCCTCTATTGATGTTACACCTGGACACCCCTTTT"
          "CTGGCATCCAGTTTTTAATCTTCAGTGGCATGTGAGATTCTCCGAAATTAATTAAAGCAAT"
          "CACACAATTCTCTCGGATACCACCTCGGTTGAAACTGACAGGTGGTTTGTTACGCATGCTA"
          "ATGCAAAGGAGCCTATATACCTTTGGCTCGGCTGCTGTAACAGGGAATATAAAGGGCAGCA"
          "TAATTTAGGAGTTTAGTGAACTTGCAACATTTACTATTTTCCCTTCTTACGTAAATATTTT"
          "TCTTTTTAATTCTAAATCAATCTTTTTCAATTTTTTGTTTGTATTCTTTTCTTGCTTAAAT"
          "CTATAACTACAAAAAACACATACATAAACTAAAA")    

    DET1 = stdgenes["DET1"]
    TPI1 = stdgenes["TPI1"]
    VMS1 = stdgenes["VMS1"]
    assert tp == str(TPI1.promoter.seq)
    assert DET1.terminator.seq == TPI1.promoter.seq
    assert TPI1.terminator.rc().seq == VMS1.terminator.seq
    assert stdgenes["TEF1"].promoter.description == 'BK006949.2 REGION: 700015..700593'
    assert stdgenes["TPI1"].promoter.description == 'BK006938.2 REGION: complement(556473..557055)'


    assert len(stdgenes["FUN26"].locus()) == 3554
    assert len(stdgenes["FUN26"].cds) == 1554
    FUN26 = stdgenes["FUN26"]
    assert FUN26.sysname == 'YAL022C'
    assert FUN26.pred.sysname == 'YAL023C'
    assert FUN26.succ.sysname == 'YAL021C'

    assert str(FUN26.promoter.seq)   in str(FUN26.locus().seq)
    assert str(FUN26.terminator.seq) in str(FUN26.locus().seq)
    assert str(stdgenes["PMT2"].promoter.seq)    in str(stdgenes["PMT2"].locus().seq)
    assert str(stdgenes["LTE1"].terminator.seq)  in str(stdgenes["LTE1"].locus().seq)

    assert str(stdgenes["DEP1"].promoter.seq) == str(stdgenes["SYN8"].promoter.seq.reverse_complement())
    assert str(stdgenes["SPO7"].promoter.seq) == str(stdgenes["MDM10"].promoter.seq.reverse_complement())

    assert str(stdgenes["FUN14"].terminator.seq) == str(stdgenes["ERP2"].terminator.seq.reverse_complement())
    assert str(stdgenes["CYS3"].promoter.seq) == str(stdgenes["DEP1"].terminator.seq)

    assert str(stdgenes["CCR4"].promoter.seq) == str(stdgenes["ATS1"].terminator.seq)


    assert stdgenes["TPI1"].tandem 
    assert not stdgenes["TPI1"].divergent 
    assert not stdgenes["GAL1"].tandem 
    assert stdgenes["GAL1"].divergent 


    assert str(stdgenes["CLN3"].promoter.seq) in str(stdgenes["CLN3"].locus(2000, 2000).seq)
    assert str(stdgenes["CLN3"].terminator.seq) in str(stdgenes["CLN3"].locus(2000, 2000).seq)

    assert str(stdgenes["CYC3"].terminator.seq) in str(stdgenes["CYC3"].locus(2000, 2000).seq)

    assert str(stdgenes["CYC3"].promoter.seq)    in str(stdgenes["CYC3"].locus(2500,2500).seq)
    assert str(stdgenes["CYC3"].terminator.seq)  == str(stdgenes["CLN3"].promoter.seq)
    assert str(stdgenes["JEN1"].promoter.seq)    == str(stdgenes["SRY1"].promoter.seq.reverse_complement())
    assert str(stdgenes["OSM1"].promoter.seq)    == str(stdgenes["ISY1"].terminator.seq)
    assert str(stdgenes["CYC1"].terminator.seq)  in str(stdgenes["CYC1"].locus().seq)
    assert str(stdgenes["UTR1"].terminator.seq)  in str(stdgenes["UTR1"].locus().seq)
    assert str(stdgenes["CDC24"].promoter.seq)   in str(stdgenes["CDC24"].locus().seq)
    assert str(stdgenes["CDC24"].terminator.seq) in str(stdgenes["CDC24"].locus().seq)

    assert stdgenes["TDH3"].gfp_cassette_kanmx.useguid() == "IJ3DVwFHTUZJq3uFO5ozwBULyME"

    s = stdgenes["CYC1"].cassette_integration_locus()

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
    
    
def test_repr():
    from unittest.mock import MagicMock
    from pygenome.saccharomyces_cerevisiae.S288C import gene_dicts
    genes, stdgenes = gene_dicts()
    s = stdgenes["CYC1"]
    pp = MagicMock()
    s._repr_pretty_(pp, None)
    pp.text.assert_any_call("Gene {}/{}".format(s.stdname, s.sysname))

    assert s._repr_html_() == "<a href='http://www.yeastgenome.org/locus/YJR048W' target='_blank'>Gene CYC1/YJR048W</a>"
    assert len(s) == 330
    assert s.short_description == "Cytochrome c, isoform 1; also known as iso-1-cytochrome c; electron carrier of mitochondrial intermembrane space that transfers electrons from ubiquinone-cytochrome c oxidoreductase to cytochrome c oxidase during cellular respiration; CYC1 has a paralog, CYC7, that arose from the whole genome duplication; human homolog CYC1 can complement yeast null mutant; mutations in human CYC1 cause insulin-responsive hyperglycemia"















    # # set local data file to be the same age remote
    # # local files should be kept
    # for fn, url in zip(_data_files,_data_urls):

    #     path = pathlib.Path( os.path.join(data_dir, fn) )
    #     os.utime(str(path), times=(path.stat().st_atime, 978307200))
    #     flo = io.BytesIO(b"some text data that will not be used")
    #     requests_mock.get(url,
    #                       headers={'last-modified'  : 'Mon, 01 Jan 2001 00:00:00 GMT', #978307200
    #                                 'content-length' : "100"},
    #                       body = flo)
    # updater()

    # # remove local data files
    # # new files should be downloaded
    # for fn, url in zip(_data_files,_data_urls):
    #     path = pathlib.Path( os.path.join(data_dir, fn) )
    #     path.unlink()
    #     flo = io.BytesIO(b"some text data that will be used")
    #     requests_mock.get(url,
    #                       headers={'last-modified'  : 'Mon, 01 Jan 2001 00:00:00 GMT', #978307200
    #                                 'content-length' : "100"},
    #                       body = flo)
    # updater()

    # # set local files newer than remote
    # # local files should be kept
    # # a critical warning should be written to log

    # for fn, url in zip(_data_files,_data_urls):
    #     path = pathlib.Path( os.path.join(data_dir, fn) )

    #     flo = io.BytesIO(b"some text data that will not be used")
    #     requests_mock.get(url,
    #                       headers={'last-modified'  : 'Sat, 01 Jan 2000 00:00:00 GMT', #978307200
    #                                 'content-length' : "100"},
    #                       body = flo)
    # updater()