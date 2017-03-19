#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''This module provides access to the Saccharomyces cerevisiae genome from Python.
   Sequences can be accessed as Bio.SeqRecord objects provided by Biopython.
'''

import re            as _re
import random        as _random
import os            as _os
import itertools     as _itertools
import urllib        as _urllib
import urllib.parse  as _urlparse
import pickle        as _pickle
import sys           as _sys
import csv           as _csv
import collections   as _collections
import time          as _time
import datetime      as _datetime

from Bio             import SeqIO           as _SeqIO
from Bio.Seq         import Seq             as _Seq
from Bio.SeqRecord   import SeqRecord       as _SeqRecord
from Bio.SeqFeature  import SeqFeature      as _SeqFeature
from Bio.SeqFeature  import FeatureLocation as _FeatureLocation

import percache as _percache
import appdirs  as _appdirs

from pydna._pretty import pretty_str as _ps
from pydna.readers import read as _read

from pygenome._pFA6a_kanMX4 import plasmid as _plasmid
_pFA6_kanMX4 = _read(_plasmid) # AJ002680

from pydna.dseqrecord import Dseqrecord as _Dseqrecord
from pydna.amplify    import pcr        as _pcr
from pydna.assembly   import Assembly   as _Assembly

def _reporthook(blocknum, blocksize, totalsize):
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

def _download(_missing_files=None):
    '''
    Download the sequence files from Saccharomyces Genome Database (www.sgd.org)
    This is typically only done once.
    '''

    #if not _missing_files:
    #    _missing_files = sorted(self._chromosome_files.values())

    _sys.stderr.write("Data files to be downloaded are:\n")

    for _file_ in _missing_files:
        _sys.stderr.write(_file_+"\n")

    _sys.stderr.write("these files will be put in {}\n".format(data_dir))

    for _file_ in sorted(_missing_files):
        _sys.stderr.write("downloading {} ".format(_file_))

        last_modified = _urllib.request.urlopen(_urlparse.urljoin(base_url, _file_)).headers['last-modified']
        remotedate = _time.strptime(last_modified, '%a, %d %b %Y %H:%M:%S %Z')

        rdate2 = int(_time.mktime(remotedate))

        _urllib.request.urlretrieve( _urlparse.urljoin(base_url, _file_),
                            filename= _os.path.join(data_dir ,_file_),
                            reporthook = _reporthook)

        _os.utime(_os.path.join(data_dir ,_file_) ,(rdate2, rdate2))

        _sys.stderr.write("{} successfully downloaded\n".format(_file_))

def _standard_name(gene):

    gene = gene.upper()

    if not _re.match("Y[A-P](R|L)\d{3}(W|C)(-.)*", gene[:7]) and gene in _feature_list:
        raise Exception("{} is not a systematic gene name.".format(gene))
    else:
        try:
            gene = _systematic_to_standard[gene]
        except KeyError:
            return None
    return gene


def _systematic_name(gene):

    gene = gene.upper()

    if _re.match("Y[A-P](R|L)\d{3}(W|C)(-.)*", gene[:7]) and gene in _feature_list:
        return _ps(gene)
    else:
        try:
            gene = _standard_to_systematic[gene]
        except KeyError:
            raise Exception("gene {} does not exist".format(gene))
    return _ps(gene)

def update():
    print("checking for updated chromosome files at\n{}".format(base_url))
    for _file_ in sorted(_chromosome_files.values()):

        last_modified = _urllib.request.urlopen(_urlparse.urljoin(base_url, _file_)).headers['last-modified']
        remotedate = _time.strptime(last_modified, '%a, %d %b %Y %H:%M:%S %Z')
        # http://stackoverflow.com/questions/5022083/how-can-i-get-the-last-modified-time-with-python3-urllib
        _missing_files = []

        if _datetime.datetime(*remotedate[:-2]) > _datetime.datetime.fromtimestamp((_os.path.getmtime(_os.path.join(data_dir, _file_)))):
            _missing_files.append(_file_)
            print("{} is available in a newer version".format(_file_))
        else:
            print("{} is the newest version".format(_file_))

        if _missing_files:
            _download(_missing_files)




if _os.getenv("CI"):
    data_dir =_os.path.join(_os.getcwd(),"DATA")
else:
    data_dir = _ps(_appdirs.user_data_dir(_os.path.join("pygenome","Saccharomyces_cerevisiae")))
    #/home/bjorn/.local/share/sgd_genome_data_files

if not _os.path.isdir(data_dir):
    _os.makedirs(data_dir)

base_url = _ps("http://downloads.yeastgenome.org/sequence/S288C_reference/NCBI_genome_source/")

_chromosome_files = { "A":"chr01.gb", "B":"chr02.gb",
                     "C":"chr03.gb", "D":"chr04.gb",
                     "E":"chr05.gb", "F":"chr06.gb",
                     "G":"chr07.gb", "H":"chr08.gb",
                     "I":"chr09.gb", "J":"chr10.gb",
                     "K":"chr11.gb", "L":"chr12.gb",
                     "M":"chr13.gb", "N":"chr14.gb",
                     "O":"chr15.gb", "P":"chr16.gb", }

_missing_files=[]

for _file_ in list(_chromosome_files.values()):
    if not _os.path.exists(_os.path.join(data_dir, _file_)):
        print("data file", _file_, "is missing")
        _missing_files.append(_file_)

if _missing_files:
    _download(_missing_files)


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
primer_url = _ps("http://www-sequence.stanford.edu/group/yeast_deletion_project/Deletion_primers_PCR_sizes.txt")
_url, _fn =_os.path.split(primer_url)

if not _os.path.exists(_os.path.join(data_dir, _fn)):
    _sys.stderr.write("\nData file {} not found\n\n".format(_fn))
    _urllib.request.urlretrieve( primer_url,
                                 _os.path.join(data_dir, _fn),
                                 reporthook = _reporthook)
    _sys.stderr.write("{} successfully downloaded and saved in {}\n".format(_fn, data_dir))

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
    _primers = _pickle.load( open(_os.path.join(data_dir, "primers.pickle"), "rb" ) )
except IOError:
    with open(_os.path.join(data_dir, _fn), 'rt') as csvfile:
        rd = _csv.reader(csvfile, delimiter='\t')
        field_names = [x.strip() for x in next(rd)]
        next(rd)
        _primers = _collections.defaultdict(tuple)
        for line_ in rd:
            v = primertuple(*[x.strip() for x in line_])
            _primers[v.ORF_name] = v
        _pickle.dump( _primers, open(_os.path.join(data_dir,"primers.pickle"), "wb" ), -1 )

_not_done_url = _ps("http://www-sequence.stanford.edu/group/yeast_deletion_project/ORFs_not_available.txt")
_url, _fn =_os.path.split(_not_done_url)
if not _os.path.exists(_os.path.join(data_dir, _fn)):
    _sys.stderr.write("\nData file {} not found\n\n".format(_fn))
    _urllib.request.urlretrieve(_not_done_url,
                       _os.path.join(data_dir, _fn),
                        reporthook = _reporthook)
    _sys.stderr.write("{} successfully downloaded and saved in {}\n".format(_fn, data_dir))

_not_done_tuple=_collections.namedtuple("not_done_tuple", "ORF_name Gene_name SGD_class")

try:
    _not_done = _pickle.load( open(_os.path.join(data_dir, "not_done.pickle"), "rb" ) )
except IOError:
    with open(_os.path.join(data_dir, _fn), 'rt') as csvfile:
        rd = _csv.reader(csvfile, delimiter='\t')
        next(rd)
        next(rd)
        next(rd)
        field_names = [x.strip() for x in next(rd)]
        next(rd)
        _not_done = _collections.defaultdict(tuple)
        for line_ in rd:
            v = _not_done_tuple(*[x.strip() for x in line_])
            _not_done[v.ORF_name] = v
        _pickle.dump( _primers, open(_os.path.join(data_dir,"not_done.pickle"), "wb" ), -1 )

try:
    _feature_list = _pickle.load( open(_os.path.join(data_dir, "feature_list.pickle"), "rb" ) )
    _standard_to_systematic = _pickle.load( open(_os.path.join(data_dir, "_standard_to_systematic.pickle"), "rb" ) )
    _systematic_to_standard = _pickle.load( open(_os.path.join(data_dir, "_systematic_to_standard.pickle"), "rb" ) )
    _systematic_to_genbank_accession = _pickle.load( open(_os.path.join(data_dir, "_systematic_to_genbank_accession.pickle"), "rb" ) )
except IOError:
    _cds          = {}
    _feature_list  = []
    _standard_to_systematic = {}
    _systematic_to_genbank_accession = {}

    for _f in list(_chromosome_files.values()):
        _krom  =  _SeqIO.read(_os.path.join(data_dir, _f), "gb")
        _features = [f for f in _krom.features if f.type=="CDS"]
        _systematic_to_genbank_accession.update( {_ps(f.qualifiers['locus_tag'][0]) : _ps("{} REGION: ".format(_krom.id)+{ 1:"{}..{}".format(f.location.start+1, f.location.end), -1:"complement({}..{})".format(f.location.start+1, f.location.end)}[f.location.strand]) for f in _features })
        _feature_list.extend( [_ps(f.qualifiers['locus_tag'][0]) for f in _features] )
        _standard_to_systematic.update( {_ps(f.qualifiers['gene'][0]):_ps(f.qualifiers['locus_tag'][0]) for f in _features if "gene" in list(f.qualifiers.keys())} )

    _systematic_to_standard = {v: k for k, v in list(_standard_to_systematic.items())}
    _pickle.dump( _feature_list, open(_os.path.join(data_dir,"feature_list.pickle"), "wb" ), -1 )
    _pickle.dump( _standard_to_systematic, open(_os.path.join(data_dir,"_standard_to_systematic.pickle"), "wb" ), -1 )
    _pickle.dump( _systematic_to_genbank_accession, open(_os.path.join(data_dir,"_systematic_to_genbank_accession.pickle"), "wb" ), -1 )
    _pickle.dump( _systematic_to_standard, open(_os.path.join(data_dir,"_systematic_to_standard.pickle"), "wb" ), -1 )

_cache = _percache.Cache(_os.path.join( data_dir, "sgd-cache" ))

del(primertuple)

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
    >>> sg.gene["TDH3"]
    YGR192C TDH3
    >>> sg.gene["TDH3"].upstream_gene
    YGR193C PDX1
    >>> sg.gene["YGR193C"].upstream_gene
    YGR194C XKS1
    >>> len(sg.intergenic_sequence("YGR192C", "YGR193C"))
    698
    >>> len(sg.intergenic_sequence("YGR192C", "YGR194C"))
    2262
    >>>

    '''

    upgene = _systematic_name(upgene)
    dngene = _systematic_name(dngene)

    if not upgene and dngene:
        raise Exception("Both upgene and dngene are needed.")
    #print upgene[1], dngene[1]
    if upgene[1] != dngene[1]:
        raise Exception("Both genes has to be on the same chromosome.")

    _krom = _SeqIO.read(_os.path.join(data_dir, _chromosome_files[upgene[1]]),"gb")
    cds  = {f.qualifiers['locus_tag'][0] :  f for f in [f for f in _krom.features if f.type=="CDS"]}
    upfeature = cds[upgene]
    startup, stopup  = upfeature.location.start,upfeature.location.end
    dnfeature = cds[dngene]
    startdn, stopdn = dnfeature.location.start, dnfeature.location.end

    assert sorted( (startup, stopup, startdn, stopdn) ) == list(_itertools.chain.from_iterable(sorted( [sorted((startup,stopup)),sorted((startdn, stopdn))] )))

    length,a,b = min([ (abs(a-b), a, b) for a, b in _itertools.product((startup,stopup),(startdn,stopdn))])

    if a<b:
        igs = _krom[a:b]
        igs.description = "{} REGION: {}..{}".format(_krom.id, a+1,b)
    else:
        igs = _krom[b:a].reverse_complement()
        igs.description = "{} REGION: complement({}..{})".format(_krom.id, b+1,a)

    igs.name = _krom.name

    igs.id = _krom.id

    igs.pydna_code = _ps("gb = pydna.Genbank('my@email.com')\n"
                         "seq = gb.nucleotide('{}')".format(igs.description))
    return igs

class _locus():

    def __init__(self, gene_name):
        self.sys = _ps(_systematic_name(gene_name))
        self.std = _ps(_standard_name(gene_name) or self.sys)
        self.chr = _ps(_chromosome_files[self.sys[1]])

    def __len__(self):
        return len(self.locus(upstream=0, downstream=0))

    def locus(self, upstream=1000, downstream=1000):

        _krom = _SeqIO.read( _os.path.join(data_dir,self.chr),"gb")

        feature = {f.qualifiers['locus_tag'][0] :  f for f in [f for f in _krom.features if f.type=="CDS"]}[self.sys]

        color = '#%02x%02x%02x' % (int(_random.uniform(150,255)),
                                   int(_random.uniform(150,255)),
                                   int(_random.uniform(150,255)),)

        feature.qualifiers.update({"ApEinfo_fwdcolor" : color,
                                   "ApEinfo_revcolor" : color,
                                   })

        start, stop = feature.location.start, feature.location.end

        if self.sys[6]=="W":
            lcs = _krom[start-upstream:stop+downstream]
            lcs.description = "{} REGION: {}..{}".format(_krom.name,start-upstream+1, stop+downstream)
        else:
            lcs = _krom[start-upstream:stop+downstream].reverse_complement()
            lcs.description = "{} REGION: complement({}..{})".format(_krom.name, stop+downstream+1, start-upstream)

        lcs.name = _krom.name

        lcs.id = _krom.id

        lcs.pydna_code = _ps("gb = pydna.Genbank('my@email.com')\n"
                             "seq = gb.nucleotide('{}')".format(lcs.description))
        return lcs

    @property
    def orf(self):
        '''
        Returns the open reading frame associated with a standard name
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
        >>> sg.gene["TDH3"].cds
        SeqRecord(seq=Seq('ATGGTTAGAGTTGCTATTAACGGTTTCGGTAGAATCGGTAGATTGGTCATGAGA...TAA', IUPACAmbiguousDNA()), id='BK006941.2', name='BK006941', description='BK006941 REGION: complement(883811..882811)', dbxrefs=[])
        >>> len(sg.gene["TDH3"].cds)
        999
        '''

        return self.locus(upstream=0, downstream=0)

    @property
    def cds(self):
        '''

        '''
        ft = [f for f in self.orf.features if f.type=="CDS"]
        ft = ft.pop()
        return ft.extract(self.orf)

    @property
    def upstream_gene(self):
        '''
        Returns the coding sequence (cds) assciated with the gene upstream
        of gene. This is defined as the gene on the chromosome located
        5' of the transcription start point of gene.
        The gene can be given as a standard name
        (eg. CYC1) or a systematic name (eg. YJR048W).

        '''

        if self.sys[6]=="W":
            sn = _feature_list[_feature_list.index(self.sys)-1]
        else:
            sn = _feature_list[_feature_list.index(self.sys)+1]
        return _locus(sn)

    @property
    def downstream_gene(self):
        '''
        Returns the coding sequence (cds) assciated with the gene downstream
        of gene. This is defined as the gene on the chromosome located
        3' of the transcription stop point of gene.
        The gene can be given as a standard name
        (eg. CYC1) or a systematic name (eg. YJR048W).

        '''

        if self.sys[6]=="C":
            sn = _feature_list[_feature_list.index(self.sys)-1]
        else:
            sn = _feature_list[_feature_list.index(self.sys)+1]
        return _locus(sn)

    @property
    def short_description(self):
        return _ps(self.locus(upstream=0, downstream=0).features[2].qualifiers["note"][0])

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
        >>> sg.gene["TDH3"].promoter
        SeqRecord(seq=Seq('ATAAAAAACACGCTTTTTCAGTTCGAGTTTATCATTATCAATACTGCCATTTCA...AAA', IUPACAmbiguousDNA()), id='YGR193C_YGR192C', name='.', description='BK006941.2 REGION: complement(883811..884508)', dbxrefs=[])
        >>> str(sg.gene["TDH3"].promoter) == str(sg.gene["YGR192C"].promoter)
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
        >>> from pygenome import sg
        >>> sg.gene["TDH3"].terminator
        SeqRecord(seq=Seq('GTGAATTTACTTTAAATCTTGCATTTAAATAAATTTTCTTTTTATAGCTTTATG...CCT', IUPACAmbiguousDNA()), id='YGR192C_YGR193C', name='.', description='Intergenic sequence between upstream gene YGR192C and downstream gene YGR191W HIP1', dbxrefs=[])
        >>> len(sg.gene["TDH3"].terminator)
        580
        >>> sg.gene["TDH3"].upstream_gene
        YGR193C PDX1
        >>> str(sg.gene["PDX1"].terminator.seq) == str(sg.gene["TDH3"].promoter.seq)
        True
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

        p = _primers[self.sys]

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

        cas = _pcr( ups, dns, _pcr( upt, dnt, _pFA6_kanMX4))

        cas.add_feature(0,len(cas), )

        locus = _Dseqrecord( self.locus(upstream, downstream) )

        asm = _Assembly( (locus, cas), max_nodes=3 )

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
    gene = _pickle.load( open( _os.path.join(data_dir, "gene.pickle"), "rb" ) )
except IOError:
    gene = {}
    for _f in _feature_list:
        gene[_f] = _locus(_f)

    for _f,_g in list(_standard_to_systematic.items()):
        gene[_f] = _locus(_g)

    _pickle.dump( gene, open( _os.path.join(data_dir,"gene.pickle"), "wb" ), -1 )

#class _g(object):
#    pass

#gn = _g()

#for _g,_l in gene.items():
#    setattr(gn,_g,_l)

#for _g,_l in _standard_to_systematic.items():
#    setattr(gn,_g,_l)


if __name__=="__main__":
    pass
    #x = gene["CYC1"]

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

