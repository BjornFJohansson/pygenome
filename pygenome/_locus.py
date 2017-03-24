#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''This module provides access to the Saccharomyces cerevisiae genome from Python.
   Sequences can be accessed as Bio.SeqRecord objects provided by Biopython.
'''

import sys           as _sys
import random        as _random
import os            as _os
import pickle        as _pickle

from Bio             import SeqIO           as _SeqIO
from Bio.Seq         import Seq             as _Seq
from Bio.SeqRecord   import SeqRecord       as _SeqRecord
from Bio.SeqFeature  import SeqFeature      as _SeqFeature
from Bio.SeqFeature  import FeatureLocation as _FeatureLocation

from pygenome._pretty import pretty_str as _ps

from pydna.readers import read as _read
from pydna.dseqrecord import Dseqrecord as _Dseqrecord
from pydna.amplify    import pcr        as _pcr
from pydna.assembly   import Assembly   as _Assembly

from pygenome._pFA6a_kanMX4 import plasmid as _plasmid
_pFA6_kanMX4 = _read(_plasmid) # AJ002680

from pygenome.systematic_name import _systematic_name
from pygenome.standard_name   import _standard_name
from pygenome._data_files     import _chromosome_files
from pygenome.intergenic      import  intergenic_sequence

data_dir = _os.path.join( _os.getenv("pygenome_data_dir"), "Saccharomyces_cerevisiae")

_feature_list = _pickle.load( open(_os.path.join(data_dir, "feature_list.pickle"), "rb" ) )
_primers = _pickle.load( open(_os.path.join(data_dir, "primers.pickle"), "rb" ) ) 


class Gene():

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
        return Gene(sn)

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
        return Gene(sn)

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
        SeqRecord(seq=Seq('GTGAATTTACTTTAAATCTTGCATTTAAATAAATTTTCTTTTTATAGCTTTATG...CCT', IUPACAmbiguousDNA()), id='YGR192C_YGR193C', name='.', description='Intergenic sequence between upstream gene YGR192C and downstream gene Gene HIP1/YGR191W', dbxrefs=[])
        >>> len(sg.gene["TDH3"].terminator)
        580
        >>> sg.gene["TDH3"].upstream_gene
        Gene PDX1/YGR193C
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
        return "Gene {}/{}".format(self.std, self.sys)

if __name__=="__main__":
    pass