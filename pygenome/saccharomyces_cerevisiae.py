#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''This module provides access to the Saccharomyces cerevisiae genome from Python.
   Sequences can be accessed as Bio.SeqRecord objects provided by Biopython.
'''
import random as _random
import os
import itertools
import urllib
import urlparse
import cPickle as pickle
import sys

from Bio             import SeqIO
from Bio.SeqFeature  import SeqFeature
from Bio.SeqFeature  import FeatureLocation

import percache
import appdirs

class pretty_str(str):
    ''' This function provide nicer output of strings in the IPython shell '''
    def _repr_pretty_(self, p, cycle):
        p.text(self)

class saccharomyces_cerevisiae_genome():

    if os.getenv("DRONE") or os.getenv("CI"):
        data_dir = os.path.join(os.getcwd(),"DATA")
    else:
        data_dir = appdirs.user_data_dir("sgd_genome_data_files")
        #/home/bjorn/.local/share/sgd_genome_data_files

    if not os.path.isdir(data_dir):
        os.mkdir(data_dir)

    cache = percache.Cache( os.path.join( data_dir, "sgd-cache" ))

    def __init__(self):

        self.data_dir = saccharomyces_cerevisiae_genome.data_dir

        self.base_url = "http://downloads.yeastgenome.org/sequence/S288C_reference/NCBI_genome_source/"

        self.chromosome_files = {   "A":"chr01.gb", "B":"chr02.gb",
                                    "C":"chr03.gb", "D":"chr04.gb",
                                    "E":"chr05.gb", "F":"chr06.gb",
                                    "G":"chr07.gb", "H":"chr08.gb",
                                    "I":"chr09.gb", "J":"chr10.gb",
                                    "K":"chr11.gb", "L":"chr12.gb",
                                    "M":"chr13.gb", "N":"chr14.gb",
                                    "O":"chr15.gb", "P":"chr16.gb",}

        missing_files=[]

        for file_ in self.chromosome_files.values():
            if not os.path.exists(os.path.join(self.data_dir, file_)):
                print "data file", file_, "is missing"
                missing_files.append(file_)

        if missing_files:
            self.download(missing_files)

        try:
            self.feature_list = pickle.load( open( os.path.join(self.data_dir, "feature_list.p"), "rb" ) )
            self.gene_to_syst = pickle.load( open( os.path.join(self.data_dir, "gene_to_syst.p"), "rb" ) )
            self.syst_to_genbank_accession = pickle.load( open( os.path.join(self.data_dir, "gene_to_genbank_accession.p"), "rb" ) )
            #a=1/0
        except:
            self._cds = {}
            self.feature_list=[]
            self.gene_to_syst={}
            #self.syst_to_syst={}
            self.syst_to_genbank_accession = {}
            for f in self.chromosome_files.values():
                krom  =  SeqIO.read(os.path.join(self.data_dir, f), "gb")
                features = [f for f in krom.features if f.type=="CDS"]  # NC_001147 REGION: complement(1071793..1072923)
                self.syst_to_genbank_accession.update( {f.qualifiers['locus_tag'][0] : "{} REGION: ".format(krom.id)+{ 1:"{}..{}".format(f.location.start+1, f.location.end), -1:"complement({}..{})".format(f.location.start+1, f.location.end)}[f.location.strand] for f in features} )
                self.feature_list.extend( [f.qualifiers['locus_tag'][0] for f in features] )
                self.gene_to_syst.update( {f.qualifiers['gene'][0]:f.qualifiers['locus_tag'][0] for f in features if "gene" in f.qualifiers.keys()} )

            pickle.dump( self.feature_list, open( os.path.join(self.data_dir,"feature_list.p"), "wb" ), -1 )
            pickle.dump( self.gene_to_syst, open( os.path.join(self.data_dir,"gene_to_syst.p"), "wb" ), -1 )
            pickle.dump( self.syst_to_genbank_accession, open( os.path.join(self.data_dir,"gene_to_genbank_accession.p"), "wb" ), -1 )

    def chromosome(self, id):
        '''
        Returns the chromosome associated with the number id

        Parameters
        ----------
        id : int or str
            chromosome number (1-16) or ("A"-"P")

        Returns
        -------
        out : Bio.SeqRecord
            Bio.SeqRecord object

        See Also
        --------
        chromosomes

        Notes
        -----
        Some of the yeast chromosomes
        return large sequences:

        ------- --------
         chr    size(bp)
        ------- --------
        chr A   230218
        chr B   316620
        chr C   813184
        chr D   576874
        chr E   1531933
        chr F   1090940
        chr G   270161
        chr H   439888
        chr I   562643
        chr J   666816
        chr K   745751
        chr L   924431
        chr M   1078177
        chr N   1091291
        chr O   784333
        ------- -------


        Examples
        --------
        >>> from pygenome import sg
        >>> len(sg.chromosome(1))
        230218
        >>> sg.chromosome(1)
        SeqRecord(seq=Seq('CCACACCACACCCACACACCCACACACCACACCACACACCACACCACACCCACA...GGG', IUPACAmbiguousDNA()), id='BK006935.2', name='BK006935', description='TPA: Saccharomyces cerevisiae S288c chromosome I.', dbxrefs=[])
        >>> len(sg.chromosome(16))
        948066
        >>> len(sg.chromosome("A"))
        230218
        >>>
        '''
        try:
            id=int(id)-1
            return SeqIO.read(os.path.join(self.data_dir, self.chromosome_files.values()[id]),"gb")
        except ValueError:
            pass
        if 1 <= (ord(id.lower())-96) <= 16:
            return SeqIO.read(os.path.join(self.data_dir, self.chromosome_files[id.upper()]),"gb")

    def chromosomes(self):
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
        return ( SeqIO.read(os.path.join(self.data_dir, f),"gb")
                 for f in self.chromosome_files.values())

    def download(self, missing_files=None):
        '''
        Download the sequence files from Saccharomyces Genome Database (www.sgd.org)
        This is typically only done once.
        '''

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

        if not missing_files:
            missing_files = sorted(self.chromosome_files.values())

        sys.stderr.write("Data files to be downloaded are:")

        for file_ in missing_files:
            sys.stderr.write(file_+"\n")

        sys.stderr.write("these files will be put in {}\n".format(self.data_dir))

        for file_ in sorted(missing_files):
            sys.stderr.write("downloading {} ".format(file_))
            urllib.urlretrieve( urlparse.urljoin(self.base_url, file_),
                                os.path.join(self.data_dir ,file_),
                                reporthook = reporthook)
            sys.stderr.write("{} successfully downloaded\n".format(file_))



    @cache
    def promoter(self, gene):
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



        gene = self.systematic_name(gene)
        if not gene:
            return
        upstream_gene = self.upstream_gene(gene)
        pr = self.intergenic_sequence(upstream_gene, gene)
        pr.features.append(SeqFeature(FeatureLocation(0, len(pr)),
                                      type = "promoter",
                                      strand = 1,
                                      qualifiers = {"note"              : "tp {} {}".format(upstream_gene,gene),
                                                    "ApEinfo_fwdcolor": "#b1e6cc",
                                                    "ApEinfo_revcolor": "#b1e681",
                                                    }))
        return pr

    def cds(self, gene):
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

        return self.locus(gene, upstream=0, downstream=0)

    @cache
    def locus(self, gene, upstream=1000, downstream=1000):
        '''
        Returns the sequence from the locus assciated with a standard name
        (eg. CYC1) or a systematic name (eg. YJR048W).

        Parameters
        ----------
        gene : str
            standard name (eg. CYC1) or a systematic name (eg. YJR048W)
        upstream : int
            Number of nucleotides before gene (default is 1000bp)
        downstream : int
            Number of nucleotides after gene (default is 1000bp)

        Returns
        -------
        out : Bio.SeqRecord
            Bio.SeqRecord object

        See Also
        --------
        cds

        Examples
        --------
        >>> from pygenome import sg
        >>> len(sg.cds("YJR048W"))
        330
        >>> len(sg.locus("YJR048W"))
        2330
        >>> len(sg.locus("YJR048W", upstream = 100, downstream = 100))
        530
        >>>

        '''
        gene = self.systematic_name(gene)
        if not gene:
            return

        krom = SeqIO.read(os.path.join(self.data_dir,self.chromosome_files[gene[1]]),"gb")

        cds ={f.qualifiers['locus_tag'][0] :  f for f in [f for f in krom.features if f.type=="CDS"]}
        feature = cds[gene]

        color = '#%02x%02x%02x' % (_random.uniform(150,255),
                                   _random.uniform(150,255),
                                   _random.uniform(150,255),)
        feature.qualifiers.update({"ApEinfo_fwdcolor" : color,
                                   "ApEinfo_revcolor" : color,
                                   })

        start, stop = feature.location.start, feature.location.end
        lcs = krom[start-upstream:stop+downstream]

        if gene[6]=="W":
            return lcs
        else:
            return lcs.reverse_complement()

    @cache
    def terminator(self, gene):
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

        gene = self.systematic_name(gene)
        if not gene:
            return
        downstream_gene = self.downstream_gene(gene)
        tm = self.intergenic_sequence(gene, downstream_gene)
        tm.features.append(SeqFeature(FeatureLocation(0,len(tm)),
                                      type = "terminator",strand = 1,
                                      qualifiers = {"note": "tp {} {}".format(gene,downstream_gene),
                                                    "ApEinfo_fwdcolor": "#b1e6cc",
                                                    "ApEinfo_revcolor": "#b1e681",
                                                    }))

        return tm

    @cache
    def intergenic_sequence(self, upgene, dngene):
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


        upgene = self.systematic_name(upgene)
        dngene = self.systematic_name(dngene)
        if not upgene and dngene and upgene[1] == dngene[1]:
            return
        krom = SeqIO.read(os.path.join(self.data_dir, self.chromosome_files[upgene[1]]),"gb")
        cds  = {f.qualifiers['locus_tag'][0] :  f for f in [f for f in krom.features if f.type=="CDS"]}
        upfeature = cds[upgene]
        startup, stopup  = upfeature.location.start,upfeature.location.end
        dnfeature = cds[dngene]
        startdn,stopdn = dnfeature.location.start, dnfeature.location.end

        assert sorted( (startup, stopup, startdn, stopdn) ) == list(itertools.chain.from_iterable(sorted( [sorted((startup,stopup)),sorted((startdn, stopdn))] )))

        length,a,b = min([ (abs(a-b), a, b) for a, b in itertools.product((startup,stopup),(startdn,stopdn))])

        if a<b:
            return krom[a:b]
        else:
            return krom[b:a].reverse_complement()

    def systematic_name(self, gene):
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
        import re
        if re.match("Y[A-P](R|L)\d{3}(W|C)(-.)*", gene[:7]) and gene in self.feature_list:
            return gene
        else:
            try:
                gene = self.gene_to_syst[gene]
            except KeyError:
                raise Exception("gene {} does not exist".format(gene))
        return gene

    def upstream_gene(self, gene):
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

        gene = self.systematic_name(gene)
        if gene[6]=="W":
            return self.feature_list[self.feature_list.index(gene)-1]
        else:
            return self.feature_list[self.feature_list.index(gene)+1]

    def downstream_gene(self, gene):
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
        gene = self.systematic_name(gene)
        if gene[6]=="C":
            return self.feature_list[self.feature_list.index(gene)-1]
        else:
            return self.feature_list[self.feature_list.index(gene)+1]

    def promoter_genbank_accession(self, gene):
        '''
        Same as the promoter_genbank method, but returns a string representing
        a portion of a Genbank file.

        Examples
        --------
        >>> from pygenome import sg
        >>> sg.promoter_genbank_accession("TDH3")
        'BK006941.2 REGION: complement(883811..884508)'
        >>>
        '''

        dngene = self.systematic_name(gene)
        if not gene:
            return
        upgene = self.systematic_name(self.upstream_gene(gene))
        if not upgene and dngene and upgene[1] == dngene[1]:
            return
        krom = SeqIO.read(os.path.join(self.data_dir, self.chromosome_files[upgene[1]]),"gb")
        cds  = {f.qualifiers['locus_tag'][0] :  f for f in [f for f in krom.features if f.type=="CDS"]}
        upfeature = cds[upgene]
        startup, stopup  = upfeature.location.start,upfeature.location.end
        dnfeature = cds[dngene]
        startdn,stopdn = dnfeature.location.start, dnfeature.location.end

        assert sorted( (startup, stopup, startdn, stopdn) ) == list(itertools.chain.from_iterable(sorted( [sorted((startup,stopup)),sorted((startdn, stopdn))] )))

        length,a,b = min([ (abs(a-b), a, b) for a, b in itertools.product((startup,stopup),(startdn,stopdn))])

        result = "{} REGION: ".format(krom.id)
        if a<b:
            return result + "{}..{}".format(a+1, b)
        else:
            return result + "complement({}..{})".format(b+1, a)

    def cds_genbank_accession(self, gene):
        '''
        Same as the cds method, but returns a string representing
        a portion of a Genbank file.

        Examples
        --------
        >>> from pygenome import sg
        >>> sg.cds_genbank_accession("TDH3")
        'BK006941.2 REGION: complement(882812..883810)'
        >>>
        '''
        return self.syst_to_genbank_accession[self.systematic_name(gene)]

    def intergenic_sequence_genbank_accession(self, upgene, dngene):
        '''
        Same as the cds method, but returns a string representing
        a portion of a Genbank file.

        Examples
        --------
        >>> from pygenome import sg
        >>> sg.intergenic_sequence_genbank_accession("YGR192C", "YGR193C")
        'BK006941.2 REGION: 883811..884508'

        '''

        upgene = self.systematic_name(upgene)
        dngene = self.systematic_name(dngene)
        if not upgene and dngene and upgene[1] == dngene[1]:
            return
        krom = SeqIO.read(os.path.join(self.data_dir, self.chromosome_files[upgene[1]]),"gb")
        cds  = {f.qualifiers['locus_tag'][0] :  f for f in [f for f in krom.features if f.type=="CDS"]}
        upfeature = cds[upgene]
        startup, stopup  = upfeature.location.start,upfeature.location.end
        dnfeature = cds[dngene]
        startdn,stopdn = dnfeature.location.start, dnfeature.location.end

        assert sorted( (startup, stopup, startdn, stopdn) ) == list(itertools.chain.from_iterable(sorted( [sorted((startup,stopup)),sorted((startdn, stopdn))] )))

        length,a,b = min([ (abs(a-b), a, b) for a, b in itertools.product((startup,stopup),(startdn,stopdn))])

        result = "{} REGION: ".format(krom.id)

        if a<b:
            return result + "{}..{}".format(a+1, b)
        else:
            return result + "complement({}..{})".format(b+1, a)

    def tandem(self, gene):
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
        return self.systematic_name(gene)[6] == self.systematic_name(self.upstream_gene(gene))[6]

    def bidirectional(self, gene):
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
        return not self.tandem(gene)

    def promoter_pydna_code(self, gene):
        return pretty_str("seq = gb.nucleotide('{}')".format(self.promoter_genbank_accession(self.systematic_name(gene))))

    def cds_pydna_code(self, gene):
        return pretty_str(("gb = pydna.Genbank('my@email.com')\n"
                           "seq = gb.nucleotide('{}')".format(self.syst_to_genbank_accession[self.systematic_name(gene)])))
    def intergenic_sequence_pydna_code(self, upgene, dngene):
        return pretty_str("seq = gb.nucleotide('{}')".format( self.intergenic_sequence_genbank_accession(upgene, dngene) ))


if __name__=="__main__":
    import doctest
    doctest.testmod()
