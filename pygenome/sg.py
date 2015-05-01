#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''This module provides access to the Saccharomyces cerevisiae genome from Python.
   Sequences can be accessed as Bio.SeqRecord objects provided by Biopython.
'''

import re
import random as _random
import os
import itertools
import urllib
import urlparse
import cPickle as pickle
import sys
import csv
import collections
import time, datetime

from Bio             import SeqIO
from Bio.Seq         import Seq
from Bio.SeqRecord   import SeqRecord
from Bio.SeqFeature  import SeqFeature
from Bio.SeqFeature  import FeatureLocation

import percache
import appdirs

import pydna

from pFA6a_kanMX4 import plasmid
pFA6_kanMX4 = pydna.read(plasmid) # AJ002680

os.environ["pydna_cache"]="refresh"

def reporthook(blocknum, blocksize, totalsize):
    readsofar = blocknum * blocksize
    if totalsize > 0:
        percent = readsofar * 1e2 / totalsize
        s = "\r%5.1f%% %*d / %d" % (
            percent, len(str(totalsize)), readsofar, totalsize)
        sys.stderr.write(s)
        if readsofar >= totalsize: # near the end
            sys.stderr.write("\n")
    else: # total size is unknown
        sys.stderr.write("read %d\n" % (readsofar,))

def download(missing_files=None):
    '''
    Download the sequence files from Saccharomyces Genome Database (www.sgd.org)
    This is typically only done once.
    '''

    #if not missing_files:
    #    missing_files = sorted(self.chromosome_files.values())

    sys.stderr.write("Data files to be downloaded are:\n")

    for file_ in missing_files:
        sys.stderr.write(file_+"\n")

    sys.stderr.write("these files will be put in {}\n".format(data_dir))

    for file_ in sorted(missing_files):
        sys.stderr.write("downloading {} ".format(file_))
        remotedate = urllib.urlopen(urlparse.urljoin(base_url, file_)).info().getdate('last-modified')


        rdate2 = int(time.mktime(remotedate))


        urllib.urlretrieve( urlparse.urljoin(base_url, file_),
                            os.path.join(data_dir ,file_),
                            reporthook = reporthook)

        os.utime(os.path.join(data_dir ,file_) ,(rdate2, rdate2))

        sys.stderr.write("{} successfully downloaded\n".format(file_))



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

    if re.match("Y[A-P](R|L)\d{3}(W|C)(-.)*", gene[:7]) and gene in feature_list:
        return gene
    else:
        try:
            gene = gene_to_syst[gene]
        except KeyError:
            raise Exception("gene {} does not exist".format(gene))
    return gene

def update():

    print "checking for updated chromosome files at {}".format(base_url)

    for file_ in chromosome_files.values():

        remotedate = urllib.urlopen(urlparse.urljoin(base_url, file_)).info().getdate('last-modified')

        missing_files = []

        if datetime.datetime(*remotedate[:-2]) > datetime.datetime.fromtimestamp((os.path.getmtime(os.path.join(data_dir, file_)))):
            missing_files.append(file_)
            print "{} is available in a newer version".format(file_)
        else:
            print "{} is the newest version".format(file_)

        if missing_files:
            download(missing_files)




if os.getenv("DRONE") or os.getenv("CI"):
    data_dir = os.path.join(os.getcwd(),"DATA")
else:
    data_dir = appdirs.user_data_dir( os.path.join("pygenome","saccharomyces_cerevisiae"))
    #/home/bjorn/.local/share/sgd_genome_data_files

if not os.path.isdir(data_dir):
    os.makedirs(data_dir)

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
    if not os.path.exists(os.path.join(data_dir, file_)):
        print "data file", file_, "is missing"
        missing_files.append(file_)

if missing_files:
    download(missing_files)

primer_url = "http://www-sequence.stanford.edu/group/yeast_deletion_project/Deletion_primers_PCR_sizes.txt"
url, fn = os.path.split(primer_url)

if not os.path.exists(os.path.join(data_dir, fn)):
    sys.stderr.write("\nData file {} not found\n\n".format(fn))
    urllib.urlretrieve( primer_url,
                        os.path.join(data_dir, fn),
                        reporthook = reporthook)
    sys.stderr.write("{} successfully downloaded and saved in {}\n".format(fn, data_dir))

primertuple = collections.namedtuple("primertuple",  '''rec_num
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
    primers = pickle.load( open( os.path.join(data_dir, "primers.p"), "rb" ) )
except IOError:
    with open(os.path.join(data_dir, fn), 'rb') as csvfile:
        rd = csv.reader(csvfile, delimiter='\t')
        field_names = [x.strip() for x in rd.next()]
        rd.next()
        primers = collections.defaultdict(tuple)
        for line_ in rd:
            v = primertuple(*[x.strip() for x in line_])
            primers[v.ORF_name] = v
        pickle.dump( primers, open( os.path.join(data_dir,"primers.p"), "wb" ), -1 )

not_done_url = "http://www-sequence.stanford.edu/group/yeast_deletion_project/ORFs_not_available.txt"
url, fn = os.path.split(not_done_url)
if not os.path.exists(os.path.join(data_dir, fn)):
    sys.stderr.write("\nData file {} not found\n\n".format(fn))
    urllib.urlretrieve( not_done_url,
                        os.path.join(data_dir, fn),
                        reporthook = reporthook)
    sys.stderr.write("{} successfully downloaded and saved in {}\n".format(fn, data_dir))

not_done_tuple=collections.namedtuple("not_done_tuple", "ORF_name Gene_name SGD_class")

try:
    not_done = pickle.load( open( os.path.join(data_dir, "not_done.p"), "rb" ) )
except IOError:
    with open(os.path.join(data_dir, fn), 'rb') as csvfile:
        rd = csv.reader(csvfile, delimiter='\t')
        rd.next()
        rd.next()
        rd.next()
        field_names = [x.strip() for x in rd.next()]
        rd.next()
        not_done = collections.defaultdict(tuple)
        for line_ in rd:
            v = not_done_tuple(*[x.strip() for x in line_])
            not_done[v.ORF_name] = v
        pickle.dump( primers, open( os.path.join(data_dir,"not_done.p"), "wb" ), -1 )

try:
    feature_list = pickle.load( open( os.path.join(data_dir, "feature_list.p"), "rb" ) )
    gene_to_syst = pickle.load( open( os.path.join(data_dir, "gene_to_syst.p"), "rb" ) )
    syst_to_genbank_accession = pickle.load( open( os.path.join(data_dir, "gene_to_genbank_accession.p"), "rb" ) )
except IOError:
    _cds = {}
    feature_list=[]
    gene_to_syst={}
    syst_to_genbank_accession = {}

    for f in chromosome_files.values():
        krom  =  SeqIO.read(os.path.join(data_dir, f), "gb")
        features = [f for f in krom.features if f.type=="CDS"]
        syst_to_genbank_accession.update( {f.qualifiers['locus_tag'][0] : "{} REGION: ".format(krom.id)+{ 1:"{}..{}".format(f.location.start+1, f.location.end), -1:"complement({}..{})".format(f.location.start+1, f.location.end)}[f.location.strand] for f in features} )
        feature_list.extend( [f.qualifiers['locus_tag'][0] for f in features] )
        gene_to_syst.update( {f.qualifiers['gene'][0]:f.qualifiers['locus_tag'][0] for f in features if "gene" in f.qualifiers.keys()} )

    pickle.dump( feature_list, open( os.path.join(data_dir,"feature_list.p"), "wb" ), -1 )
    pickle.dump( gene_to_syst, open( os.path.join(data_dir,"gene_to_syst.p"), "wb" ), -1 )
    pickle.dump( syst_to_genbank_accession, open( os.path.join(data_dir,"gene_to_genbank_accession.p"), "wb" ), -1 )

cache = percache.Cache( os.path.join( data_dir, "sgd-cache" ))

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
    return ( SeqIO.read(os.path.join(self.data_dir, f),"gb") for f in self.chromosome_files.values())

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
    krom = SeqIO.read(os.path.join(data_dir, chromosome_files[upgene[1]]),"gb")
    cds  = {f.qualifiers['locus_tag'][0] :  f for f in [f for f in krom.features if f.type=="CDS"]}
    upfeature = cds[upgene]
    startup, stopup  = upfeature.location.start,upfeature.location.end
    dnfeature = cds[dngene]
    startdn, stopdn = dnfeature.location.start, dnfeature.location.end

    assert sorted( (startup, stopup, startdn, stopdn) ) == list(itertools.chain.from_iterable(sorted( [sorted((startup,stopup)),sorted((startdn, stopdn))] )))

    length,a,b = min([ (abs(a-b), a, b) for a, b in itertools.product((startup,stopup),(startdn,stopdn))])

    if a<b:
        igs = krom[a:b]
        igs.description = "{} REGION: {}..{}".format(krom.name, a+1,b)
    else:
        igs = krom[b:a].reverse_complement()
        igs.description = "{} REGION: complement({}..{})".format(krom.name, b+1,a)


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
        self.n   = gene_name
        self.sys = systematic_name(gene_name)
        self.chr = chromosome_files[self.sys[1]]

    def locus(self, upstream=1000, downstream=1000):

        krom = SeqIO.read( os.path.join(data_dir,self.chr),"gb")

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

        pr = intergenic_sequence(self.upstream_gene().sys, self.sys)
        pr.id = self.upstream_gene().sys+"_"+self.sys
        pr.name = "."
        #pr.description = "Intergenic sequence between upstream gene {} and downstream gene {}".format(self.upstream_gene().sys, self.sys)
        pr.features.append(SeqFeature(FeatureLocation(0, len(pr)),
                                      type = "promoter",
                                      strand = 1,
                                      qualifiers = {"note"              : "tp {} {}".format(self.upstream_gene().sys, self.sys),
                                                    "ApEinfo_fwdcolor": "#b1e6cc",
                                                    "ApEinfo_revcolor": "#b1e681",
                                                    }))
        return pr

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


        tm = intergenic_sequence( self.sys, self.downstream_gene().sys )
        tm.id = self.sys + "_" + self.upstream_gene().sys
        tm.name = "."
        tm.description = "Intergenic sequence between upstream gene {} and downstream gene {}".format( self.sys, self.downstream_gene() )
        tm.features.append(SeqFeature(FeatureLocation(0, len(tm)),
                                      type = "terminator",
                                      strand = 1,
                                      qualifiers = {"note"              : "tp {} {}".format( self.sys, self.downstream_gene() ),
                                                    "ApEinfo_fwdcolor": "#b1e6cc",
                                                    "ApEinfo_revcolor": "#b1e681",
                                                    }))
        return tm

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
        return self.sys[6] == self.upstream_gene().sys[6]

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
        return not self.tandem()

    def deletion_locus(self, upstream=1000, downstream=1000):

        p = primers[self.sys]

        if not p:
            return "No deletion primers available!"

        upt = SeqRecord( Seq( p.UPTAG_primer_sequence ))
        dnt = SeqRecord( Seq( p.DNTAG_primer_sequence ))
        ups = SeqRecord( Seq( p.UPstream45_primer_sequence ))
        dns = SeqRecord( Seq( p.DNstream45_primer_sequence ))

        if "" in [str(x.seq).strip() for x in [upt,dnt,ups,dns]]:
            return "One deletion primer missing!"

        upt.id  = "UPTAG_primer_{}".format(self.sys)
        dnt.id  = "DNTAG_primer_{}".format(self.sys)
        ups.id  = "UPstream45_{}".format(self.sys)
        dns.id  = "DNstream45_{}".format(self.sys)

        cas = pydna.pcr( ups, dns, pydna.pcr( upt, dnt, pFA6_kanMX4))

        cas.add_feature(0,len(cas), )

        locus = pydna.Dseqrecord( self.locus(upstream, downstream) )

        asm = pydna.Assembly( (locus, cas), max_nodes=3 )

        candidate = asm.linear_products[0]

        #print candidate.figure()

        kanmx4_gene = candidate

        kanmx4_gene.name = "{}::KanMX4".format(self.sys.lower())

        kanmx4_gene.id = "{} locus with {} bp up and {} bp downstream DNA".format(kanmx4_gene.name, upstream, downstream)

        #kanmx4_gene.version = "."

        k = SeqRecord( Seq( str(kanmx4_gene.seq), kanmx4_gene.seq.alphabet),
                      id = kanmx4_gene.id,
                      name = kanmx4_gene.name,
                      description = kanmx4_gene.description,
                      dbxrefs = kanmx4_gene.dbxrefs,
                      features = kanmx4_gene.features,
                      annotations = kanmx4_gene.annotations)

        return k


    def __repr__(self):
        return "yeast gene "+self.sys

try:
    gene = pickle.load( open( os.path.join(data_dir, "gene.p"), "rb" ) )
except IOError:
    gene = {}
    for f in feature_list:
        gene[f] = _locus(f)

    for f,g in gene_to_syst.items():
        gene[f] = _locus(g)

    pickle.dump( gene, open( os.path.join(data_dir,"gene.p"), "wb" ), -1 )



if __name__=="__main__":

    update()

    #import doctest
    #doctest.testmod()
    #hej = gene["YGR194C"]

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

