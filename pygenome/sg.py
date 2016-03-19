#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''This module provides access to the Saccharomyces cerevisiae genome from Python.
   Sequences can be accessed as Bio.SeqRecord objects provided by Biopython.
'''

import re as _re
import random as _random
import os as _os
import itertools as _itertools
import urllib as _urllib
import urlparse as _urlparse
import cPickle as _pickle
import sys as _sys
import csv as _csv
import collections as _collections
import time as _time
import datetime as _datetime

from Bio             import SeqIO      as _SeqIO
from Bio.Seq         import Seq        as _Seq
from Bio.SeqRecord   import SeqRecord  as _SeqRecord
from Bio.SeqFeature  import SeqFeature as _SeqFeature
from Bio.SeqFeature  import FeatureLocation as _FeatureLocation

import percache as _percache
import appdirs  as _appdirs

import pydna as _pydna

from _pFA6a_kanMX4 import plasmid as _plasmid
_pFA6_kanMX4 = _pydna.read(_plasmid) # AJ002680

#os.environ["pydna_cache"]="refresh"


def reporthook(blocknum, blocksize, totalsize):
    readsofar = blocknum * blocksize
    if totalsize > 0:
        percent = readsofar * 1e2 / totalsize
        s = "\r%5.1f%% %*d / %d" % (
            percent, len(str(totalsize)), readsofar, totalsize)
        _sys.stderr.write(s)
        if readsofar >= totalsize: # near the end
            _sys.stderr.write("\n")
    else: # total size is unknown
        _sys.stderr.write("read %d\n" % (readsofar,))

def download(missing_files=None):
    '''
    Download the sequence files from Saccharomyces Genome Database (www.sgd.org)
    This is typically only done once.
    '''

    #if not missing_files:
    #    missing_files = sorted(self.chromosome_files.values())

    _sys.stderr.write("Data files to be downloaded are:\n")

    for file_ in missing_files:
        _sys.stderr.write(file_+"\n")

    _sys.stderr.write("these files will be put in {}\n".format(data_dir))

    for file_ in sorted(missing_files):
        _sys.stderr.write("downloading {} ".format(file_))
        remotedate = _urllib.urlopen(_urlparse.urljoin(base_url, file_)).info().getdate('last-modified')


        rdate2 = int(_time.mktime(remotedate))


        _urllib.urlretrieve( _urlparse.urljoin(base_url, file_),
                            _os.path.join(data_dir ,file_),
                            reporthook = reporthook)

        _os.utime(os.path.join(data_dir ,file_) ,(rdate2, rdate2))

        _sys.stderr.write("{} successfully downloaded\n".format(file_))

def standard_name(gene):
    '''
    Returns the standard name associated with a systematic name (if any).

    Parameters
    ----------
    gene : str
        systematic name (eg. YBR020W or YGR192C)

    Returns
    -------
    out : str
        String

    Examples
    --------
    >>> from pygenome import sg
    >>> sg.standard_name("YBR020W")
    'GAL1'
    >>> sg.standard_name("YJR048W")
    'CYC1'
    >>> sg.standard_name("YGR192C")
    'TDH3'
    >>>
    '''
    gene = gene.upper()

    if not _re.match("Y[A-P](R|L)\d{3}(W|C)(-.)*", gene[:7]) and gene in feature_list:
        raise Exception("{} is not a systematic gene name.".format(gene))
    else:
        try:
            gene = syst_to_stand[gene]
        except KeyError:
            return None
    return gene


def systematic_name(gene):
    '''
    Returns the systematic name associated with a standard name.

    Parameters
    ----------
    gene : str
        standard name (eg. CYC1 or TDH3)

    Returns
    -------
    out : str
        String

    Examples
    --------
    >>> from pygenome import sg
    >>> sg.systematic_name("GAL1")
    'YBR020W'
    >>> sg.systematic_name("CYC1")
    'YJR048W'
    >>> sg.systematic_name("TDH3")
    'YGR192C'
    >>>
    '''

    gene = gene.upper()

    if _re.match("Y[A-P](R|L)\d{3}(W|C)(-.)*", gene[:7]) and gene in feature_list:
        return gene
    else:
        try:
            gene = stand_to_syst[gene]
        except KeyError:
            raise Exception("gene {} does not exist".format(gene))
    return gene

def update():
    print "checking for updated chromosome files at {}".format(base_url)
    for file_ in chromosome_files.values():

        remotedate = _urllib.urlopen(_urlparse.urljoin(base_url, file_)).info().getdate('last-modified')

        missing_files = []

        if _datetime.datetime(*remotedate[:-2]) > _datetime.datetime.fromtimestamp((os.path.getmtime(os.path.join(data_dir, file_)))):
            missing_files.append(file_)
            print "{} is available in a newer version".format(file_)
        else:
            print "{} is the newest version".format(file_)

        if missing_files:
            download(missing_files)




if _os.getenv("DRONE") or _os.getenv("CI"):
    data_dir =_os.path.join(os.getcwd(),"DATA")
else:
    data_dir = _appdirs.user_data_dir(_os.path.join("pygenome","saccharomyces_cerevisiae"))
    #/home/bjorn/.local/share/sgd_genome_data_files

if not _os.path.isdir(data_dir):
    _os.makedirs(data_dir)

base_url = "http://downloads.yeastgenome.org/sequence/S288C_reference/NCBI_genome_source/"

chromosome_files = { "A":"chr01.gb", "B":"chr02.gb",
                     "C":"chr03.gb", "D":"chr04.gb",
                     "E":"chr05.gb", "F":"chr06.gb",
                     "G":"chr07.gb", "H":"chr08.gb",
                     "I":"chr09.gb", "J":"chr10.gb",
                     "K":"chr11.gb", "L":"chr12.gb",
                     "M":"chr13.gb", "N":"chr14.gb",
                     "O":"chr15.gb", "P":"chr16.gb", }

missing_files=[]

for file_ in chromosome_files.values():
    if not _os.path.exists(_os.path.join(data_dir, file_)):
        print "data file", file_, "is missing"
        missing_files.append(file_)

if missing_files:
    download(missing_files)


'''
gn.JEN1.sgdpage
       .shortdescription
       .longdescription

Standard Name
GDH3
Systematic Name
YAL062W
SGD ID
S000000058
Aliases
FUN51
Feature Type
ORF , Verified
Description
NADP(+)-dependent glutamate dehydrogenase; synthesizes glutamate from ammonia and alpha-ketoglutarate; rate of alpha-ketoglutarate utilization differs from Gdh1p; expression regulated by nitrogen and carbon sources; GDH3 has a paralog, GDH1, that arose from the whole genome duplication 1 2 3 4
Name Description
Glutamate DeHydrogenase
Paralog
GDH1 1


'''
primer_url = "http://www-sequence.stanford.edu/group/yeast_deletion_project/Deletion_primers_PCR_sizes.txt"
url, fn =_os.path.split(primer_url)

if not _os.path.exists(_os.path.join(data_dir, fn)):
    _sys.stderr.write("\nData file {} not found\n\n".format(fn))
    _urllib.urlretrieve( primer_url,
                       _os.path.join(data_dir, fn),
                        reporthook = reporthook)
    _sys.stderr.write("{} successfully downloaded and saved in {}\n".format(fn, data_dir))

primertuple = _collections.namedtuple("primertuple", '''rec_num
                                                        ORF_name
                                                        deletion_alias
                                                        essential
                                                        A_confirmation_primer_sequence
                                                        B_confirmation_primer_sequence
                                                        C_confirmation_primer_sequence
                                                        D_confirmation_primer_sequence
                                                        UPTAG_primer_sequence
                                                        DNTAG_primer_sequence
                                                        UPstream45_primer_sequence
                                                        DNstream45_primer_sequence
                                                        UPstream90_primer_sequence
                                                        DNstream90_primer_sequence
                                                        AB_wt_PCR
                                                        AkanB_del_PCR
                                                        CD_wt_PCR
                                                        DkanC_del_PCR
                                                        AD_wt
                                                        AD_del
                                                        AB_del_PCR
                                                        CD_del_PCR
                                                        UPTAG_sequence_20mer
                                                        DNTAG_sequence_20mer''')



try:
    primers = _pickle.load( open(_os.path.join(data_dir, "primers.p"), "rb" ) )
except IOError:
    with open(_os.path.join(data_dir, fn), 'rb') as csvfile:
        rd = _csv.reader(csvfile, delimiter='\t')
        field_names = [x.strip() for x in rd.next()]
        rd.next()
        primers = _collections.defaultdict(tuple)
        for line_ in rd:
            v = primertuple(*[x.strip() for x in line_])
            primers[v.ORF_name] = v
        _pickle.dump( primers, open(_os.path.join(data_dir,"primers.p"), "wb" ), -1 )

not_done_url = "http://www-sequence.stanford.edu/group/yeast_deletion_project/ORFs_not_available.txt"
url, fn =_os.path.split(not_done_url)
if not _os.path.exists(_os.path.join(data_dir, fn)):
    _sys.stderr.write("\nData file {} not found\n\n".format(fn))
    _urllib.urlretrieve( not_done_url,
                       _os.path.join(data_dir, fn),
                        reporthook = reporthook)
    _sys.stderr.write("{} successfully downloaded and saved in {}\n".format(fn, data_dir))

not_done_tuple=_collections.namedtuple("not_done_tuple", "ORF_name Gene_name SGD_class")

try:
    not_done = _pickle.load( open(_os.path.join(data_dir, "not_done.p"), "rb" ) )
except IOError:
    with open(os.path.join(data_dir, fn), 'rb') as csvfile:
        rd = _csv.reader(csvfile, delimiter='\t')
        rd.next()
        rd.next()
        rd.next()
        field_names = [x.strip() for x in rd.next()]
        rd.next()
        not_done = _collections.defaultdict(tuple)
        for line_ in rd:
            v = not_done_tuple(*[x.strip() for x in line_])
            not_done[v.ORF_name] = v
        _pickle.dump( primers, open(_os.path.join(data_dir,"not_done.p"), "wb" ), -1 )

try:
    feature_list = _pickle.load( open(_os.path.join(data_dir, "feature_list.p"), "rb" ) )
    stand_to_syst = _pickle.load( open(_os.path.join(data_dir, "stand_to_syst.p"), "rb" ) )
    syst_to_stand = _pickle.load( open(_os.path.join(data_dir, "syst_to_stand.p"), "rb" ) )
    syst_to_genbank_accession = _pickle.load( open(_os.path.join(data_dir, "syst_to_genbank_accession.p"), "rb" ) )
except IOError:
    _cds          = {}
    feature_list  = []
    stand_to_syst = {}
    syst_to_genbank_accession = {}

    for f in chromosome_files.values():
        krom  =  _SeqIO.read(_os.path.join(data_dir, f), "gb")
        features = [f for f in krom.features if f.type=="CDS"]
        syst_to_genbank_accession.update( {f.qualifiers['locus_tag'][0] : "{} REGION: ".format(krom.id)+{ 1:"{}..{}".format(f.location.start+1, f.location.end), -1:"complement({}..{})".format(f.location.start+1, f.location.end)}[f.location.strand] for f in features} )
        feature_list.extend( [f.qualifiers['locus_tag'][0] for f in features] )
        stand_to_syst.update( {f.qualifiers['gene'][0]:f.qualifiers['locus_tag'][0] for f in features if "gene" in f.qualifiers.keys()} )

    syst_to_stand = {v: k for k, v in stand_to_syst.items()}
    _pickle.dump( feature_list, open(_os.path.join(data_dir,"feature_list.p"), "wb" ), -1 )
    _pickle.dump( stand_to_syst, open(_os.path.join(data_dir,"stand_to_syst.p"), "wb" ), -1 )
    _pickle.dump( syst_to_genbank_accession, open(_os.path.join(data_dir,"syst_to_genbank_accession.p"), "wb" ), -1 )
    _pickle.dump( syst_to_stand, open(_os.path.join(data_dir,"syst_to_stand.p"), "wb" ), -1 )

cache = _percache.Cache(_os.path.join( data_dir, "sgd-cache" ))

def chromosomes():
    ''' Returns a generator containing all yeast chromosomes in the form of Bio.SeqRecord objects

    Returns
    -------
    out : Generator
        Generator of Bio.SeqRecord object

    See Also
    --------
    chromosome

    Examples
    --------
    >>> from pygenome import sg
    >>> sg.chromosomes().next()
    SeqRecord(seq=Seq('CCACACCACACCCACACACCCACACACCACACCACACACCACACCACACCCACA...GGG', IUPACAmbiguousDNA()), id='BK006935.2', name='BK006935', description='TPA: Saccharomyces cerevisiae S288c chromosome I.', dbxrefs=[])
    >>>
    '''
    return ( _SeqIO.read(os.path.join(self.data_dir, f),"gb") for f in self.chromosome_files.values())

def intergenic_sequence(upgene, dngene):
    '''
    Returns the intergenic sequence between two genes on the same chromosome.

    Parameters
    ----------
    upgene : str
        standard name (eg. CYC1) or a systematic name (eg. YJR048W)
    dngene : str
        standard name (eg. CYC1) or a systematic name (eg. YJR048W)

    Returns
    -------
    out : Bio.SeqRecord
        Bio.SeqRecord object

    See Also
    --------
    intergenic_sequence_genbank_accession

    Examples
    --------
    >>> from pygenome import sg
    >>> sg.systematic_name("TDH3")
    'YGR192C'
    >>> sg.upstream_gene("TDH3")
    'YGR193C'
    >>> sg.upstream_gene("YGR193C")
    'YGR194C'
    >>> len(sg.intergenic_sequence("YGR192C", "YGR193C"))
    698
    >>> len(sg.intergenic_sequence("YGR192C", "YGR194C"))
    2262
    >>>

    '''

    upgene = systematic_name(upgene)
    dngene = systematic_name(dngene)

    if not upgene and dngene and upgene[1] == dngene[1]:
        return
    krom = _SeqIO.read(_os.path.join(data_dir, chromosome_files[upgene[1]]),"gb")
    cds  = {f.qualifiers['locus_tag'][0] :  f for f in [f for f in krom.features if f.type=="CDS"]}
    upfeature = cds[upgene]
    startup, stopup  = upfeature.location.start,upfeature.location.end
    dnfeature = cds[dngene]
    startdn, stopdn = dnfeature.location.start, dnfeature.location.end

    assert sorted( (startup, stopup, startdn, stopdn) ) == list(_itertools.chain.from_iterable(sorted( [sorted((startup,stopup)),sorted((startdn, stopdn))] )))

    length,a,b = min([ (abs(a-b), a, b) for a, b in _itertools.product((startup,stopup),(startdn,stopdn))])

    if a<b:
        igs = krom[a:b]
        igs.description = "{} REGION: {}..{}".format(krom.id, a+1,b)
    else:
        igs = krom[b:a].reverse_complement()
        igs.description = "{} REGION: complement({}..{})".format(krom.id, b+1,a)

    igs.name = krom.name

    igs.id = krom.id

    igs.pydna_code = pretty_str("gb = pydna.Genbank('my@email.com')\n"
                                "seq = gb.nucleotide('{}')".format(igs.description))
    return igs

class pretty_str(str):
    ''' This function provide nicer output of strings in the IPython shell '''
    def _repr_pretty_(self, p, cycle):
        p.text(self)

class _locus():

    def __init__(self, gene_name):
        self.sys = systematic_name(gene_name)
        self.std = standard_name(gene_name) or self.sys
        self.chr = chromosome_files[self.sys[1]]

    def __len__(self):
        return len(self.locus(upstream=0, downstream=0))

    def locus(self, upstream=1000, downstream=1000):

        krom = _SeqIO.read( _os.path.join(data_dir,self.chr),"gb")

        cds = {f.qualifiers['locus_tag'][0] :  f for f in [f for f in krom.features if f.type=="CDS"]}

        feature = cds[self.sys]

        color = '#%02x%02x%02x' % (_random.uniform(150,255),
                                   _random.uniform(150,255),
                                   _random.uniform(150,255),)
        feature.qualifiers.update({"ApEinfo_fwdcolor" : color,
                                   "ApEinfo_revcolor" : color,
                                   })

        start, stop = feature.location.start, feature.location.end

        if self.sys[6]=="W":
            lcs = krom[start-upstream:stop+downstream]
            lcs.description = "{} REGION: {}..{}".format(krom.name,start-upstream+1, stop+downstream)
        else:
            lcs = krom[start-upstream:stop+downstream].reverse_complement()
            lcs.description = "{} REGION: complement({}..{})".format(krom.name, stop+downstream+1, start-upstream)

        lcs.name = krom.name


        lcs.id = krom.id

        lcs.pydna_code = pretty_str("gb = pydna.Genbank('my@email.com')\n"
                                    "seq = gb.nucleotide('{}')".format(lcs.description))
        return lcs

    @property
    def cds(self):
        '''
        Returns the coding sequence assciated with a standard name
        (eg. CYC1) or a systematic name (eg. YJR048W).

        Parameters
        ----------
        gene : str
            standard name (eg. CYC1) or a systematic name (eg. YJR048W)

        Returns
        -------
        out : Bio.SeqRecord
            Bio.SeqRecord object

        See Also
        --------
        cds_genbank_accession
        cds_pydna_code

        Examples
        --------
        >>> from pygenome import sg
        >>> sg.cds("TDH3")
        SeqRecord(seq=Seq('ATGGTTAGAGTTGCTATTAACGGTTTCGGTAGAATCGGTAGATTGGTCATGAGA...TAA', IUPACAmbiguousDNA()), id='<unknown id>', name='<unknown name>', description='<unknown description>', dbxrefs=[])
        >>> len(sg.cds("TDH3"))
        999
        >>> sg.cds("YJR048W")
        SeqRecord(seq=Seq('ATGACTGAATTCAAGGCCGGTTCTGCTAAGAAAGGTGCTACACTTTTCAAGACT...TAA', IUPACAmbiguousDNA()), id='BK006943.2', name='BK006943', description='TPA: Saccharomyces cerevisiae S288c chromosome X.', dbxrefs=[])
        >>> len(sg.cds("YJR048W"))
        330
        >>>

        '''
        return self.locus(upstream=0, downstream=0)


    @property
    def upstream_gene(self):
        '''
        Returns the coding sequence (cds) assciated with the gene upstream
        of gene. This is defined as the gene on the chromosome located
        5' of the transcription start point of gene.
        The gene can be given as a standard name
        (eg. CYC1) or a systematic name (eg. YJR048W).

        Examples
        --------
        >>> from pygenome import sg
        >>> sg.systematic_name("RFA1")
        'YAR007C'
        >>> sg.upstream_gene("RFA1")
        'YAR008W'
        >>> sg.systematic_name("CYC3")
        'YAL039C'
        >>> sg.systematic_name("CYC3")
        'YAL039C'
        >>>
        '''

        if self.sys[6]=="W":
            sn = feature_list[feature_list.index(self.sys)-1]
        else:
            sn = feature_list[feature_list.index(self.sys)+1]
        return _locus(sn)

    @property
    def downstream_gene(self):
        '''
        Returns the coding sequence (cds) assciated with the gene downstream
        of gene. This is defined as the gene on the chromosome located
        3' of the transcription stop point of gene.
        The gene can be given as a standard name
        (eg. CYC1) or a systematic name (eg. YJR048W).

        Examples
        --------
        >>> from pygenome import sg
        >>> sg.downstream_gene("RFA1")
        'YAR003W'
        >>> sg.systematic_name("RFA1")
        'YAR007C'
        >>> sg.downstream_gene("CYC3")
        'YAL040C'
        >>> sg.systematic_name("CYC3")
        'YAL039C'
        >>>
        '''

        if self.sys[6]=="C":
            sn = feature_list[feature_list.index(self.sys)-1]
        else:
            sn = feature_list[feature_list.index(self.sys)+1]
        return _locus(sn)

    @property
    def short_description(self):
        return self.locus(upstream=0, downstream=0).features[2].qualifiers["note"][0]

    @property
    def promoter(self):
        '''
        Returns the sequence of the promoter assciated with
        a standard name (eg. CYC1) or a systematic name (eg. YJR048W).

        The promoter is defined as the sequence between the start codon
        of the gene and the nearest start or stop codon of the upstream
        gene.

        Parameters
        ----------
        gene : str
            standard name (eg. CYC1) or a systematic name (eg. YJR048W)

        Returns
        -------
        out : Bio.SeqRecord
            Bio.SeqRecord object

        See Also
        --------
        promoter_genbank_accession
        promoter_pydna_code


        Examples
        --------
        >>> from pygenome import sg
        >>> sg.promoter("TDH3")
        SeqRecord(seq=Seq('ATAAAAAACACGCTTTTTCAGTTCGAGTTTATCATTATCAATACTGCCATTTCA...AAA', IUPACAmbiguousDNA()), id='<unknown id>', name='<unknown name>', description='<unknown description>', dbxrefs=[])
        >>> str(sg.promoter("TDH3")) == str(sg.promoter("YGR192C"))
        True
        >>>

        '''

        pr = intergenic_sequence(self.upstream_gene.sys, self.sys)
        pr.id = self.upstream_gene.sys+"_"+self.sys
        pr.name = "."
        #pr.description = "Intergenic sequence between upstream gene {} and downstream gene {}".format(self.upstream_gene().sys, self.sys)
        pr.features.append(_SeqFeature(_FeatureLocation(0, len(pr)),
                                      type = "promoter",
                                      strand = 1,
                                      qualifiers = {"note"              : "tp {} {}".format(self.upstream_gene.sys, self.sys),
                                                    "ApEinfo_fwdcolor": "#b1e6cc",
                                                    "ApEinfo_revcolor": "#b1e681",
                                                    }))
        return pr



    @property
    def terminator(self):
        '''
        Returns the sequence of the terminator assciated with
        a standard name (eg. CYC1) or a systematic name (eg. YJR048W).

        The promoter is defined as the sequence between the stop codon
        of the gene and the nearest start or stop codon of the downstream
        gene.

        Parameters
        ----------
        gene : str
            standard name (eg. CYC1) or a systematic name (eg. YJR048W)

        Returns
        -------
        out : Bio.SeqRecord
            Bio.SeqRecord object

        See Also
        --------
        terminator_genbank_accession
        terinator_pydna_code


        Examples
        --------
        >>> sg.terminator("TDH3")
        SeqRecord(seq=Seq('GTGAATTTACTTTAAATCTTGCATTTAAATAAATTTTCTTTTTATAGCTTTATG...CCT', IUPACAmbiguousDNA()), id='<unknown id>', name='<unknown name>', description='<unknown description>', dbxrefs=[])
        >>> len(sg.terminator("TDH3"))
        580
        >>> sg.upstream_gene("TDH3")
        'YGR193C'
        >>> str(sg.terminator("YGR193C")) == str(sg.promoter("TDH3"))

        '''


        tm = intergenic_sequence( self.sys, self.downstream_gene.sys )
        tm.id = self.sys + "_" + self.upstream_gene.sys
        tm.name = "."
        tm.description = "Intergenic sequence between upstream gene {} and downstream gene {}".format( self.sys, self.downstream_gene)
        tm.features.append(_SeqFeature(_FeatureLocation(0, len(tm)),
                                      type = "terminator",
                                      strand = 1,
                                      qualifiers = {"note"              : "tp {} {}".format( self.sys, self.downstream_gene ),
                                                    "ApEinfo_fwdcolor": "#b1e6cc",
                                                    "ApEinfo_revcolor": "#b1e681",
                                                    }))

        return tm

    @property
    def tandem(self):
        '''
        Returns True if a gene is expressed in the same direction on the chromosome as
        the gene immediatelly upstream.

        ::
            Tandem genes:

            Gene1 Gene2
            ----> ---->

            bidirectional genes:

            Gene1 Gene2
            ----> <----


        Parameters
        ----------
        gene : str
            standard name (eg. CYC1 or TDH3)

        Returns
        -------
        out : bool
            Boolean; True or False

        '''
        return self.sys[6] == self.upstream_gene.sys[6]


    @property
    def bidirectional(self):
        '''
        Returns True if a gene is not expressed in the same direction on the chromosome as
        the gene immediatelly upstream.

        ::
            Tandem genes:

            Gene1 Gene2
            ----> ---->

            bidirectional genes:

            Gene1 Gene2
            ----> <----

            Gene1 Gene2
            <---- ---->


        Parameters
        ----------
        gene : str
            standard name (eg. CYC1 or TDH3)

        Returns
        -------
        out : bool
            Boolean; True or False

        '''
        return not self.tandem

    @property
    def deletion_locus(self, upstream=1000, downstream=1000):

        p = primers[self.sys]

        if not p:
            return "No deletion primers available!"

        upt = _SeqRecord( _Seq( p.UPTAG_primer_sequence ))
        dnt = _SeqRecord( _Seq( p.DNTAG_primer_sequence ))
        ups = _SeqRecord( _Seq( p.UPstream45_primer_sequence ))
        dns = _SeqRecord( _Seq( p.DNstream45_primer_sequence ))

        if "" in [str(x.seq).strip() for x in [upt,dnt,ups,dns]]:
            return "One deletion primer missing!"

        upt.id  = "UPTAG_primer_{}".format(self.sys)
        dnt.id  = "DNTAG_primer_{}".format(self.sys)
        ups.id  = "UPstream45_{}".format(self.sys)
        dns.id  = "DNstream45_{}".format(self.sys)

        cas = _pydna.pcr( ups, dns, _pydna.pcr( upt, dnt, _pFA6_kanMX4))

        cas.add_feature(0,len(cas), )

        locus = _pydna.Dseqrecord( self.locus(upstream, downstream) )

        asm = _pydna.Assembly( (locus, cas), max_nodes=3 )

        candidate = asm.linear_products[0]

        #print candidate.figure()

        kanmx4_gene = candidate

        kanmx4_gene.name = "{}::KanMX4".format(self.sys.lower())

        kanmx4_gene.id = "{} locus with {} bp up and {} bp downstream DNA".format(kanmx4_gene.name, upstream, downstream)

        #kanmx4_gene.version = "."

        k = _SeqRecord( _Seq( str(kanmx4_gene.seq), kanmx4_gene.seq.alphabet),
                      id = kanmx4_gene.id,
                      name = kanmx4_gene.name,
                      description = kanmx4_gene.description,
                      dbxrefs = kanmx4_gene.dbxrefs,
                      features = kanmx4_gene.features,
                      annotations = kanmx4_gene.annotations)

        return k

    def __repr__(self):
        return "{} {}".format(self.sys, self.std)

try:
    gene = _pickle.load( open( _os.path.join(data_dir, "gene.p"), "rb" ) )
except IOError:
    gene = {}
    for _f in feature_list:
        gene[_f] = _locus(_f)

    for _f,_g in stand_to_syst.items():
        gene[_f] = _locus(_g)

    _pickle.dump( gene, open( _os.path.join(data_dir,"gene.p"), "wb" ), -1 )



class _g(object):
    pass

gn = _g()

for _g,_l in gene.items():
    setattr(gn,_g,_l)

for _g,_l in stand_to_syst.items():
    setattr(gn,_g,_l)


if __name__=="__main__":

    x = gene["CYC1"]

    #update()

    #import doctest
    #doctest.testmod()
    #hej = gene["YGR194C"]

    #print hej.terminator()

    #print hej.deletion_locus().format("gb")

    #from pydna_helper import ape

    #ape( hej.deletion_locus() )

    #print hej.promoter().pydna_code


    #print hej.locus().format("gb")

    #print hej.locus().accession

    #print len(hej.locus())

    #print hej.upstream_gene()

    #print hej.downstream_gene()

    #print intergenic_sequence(hej.upstream_gene(), hej)

    #print hej.terminator().pydna_code

    #genes = {}

