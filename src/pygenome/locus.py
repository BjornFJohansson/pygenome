#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''This module provides access to the Saccharomyces cerevisiae genome from Python.
   Sequences can be accessed as Bio.SeqRecord objects provided by Biopython.
'''

import sys as _sys
from warnings import warn
import random        as _random
import os            as _os
import collections   as  _collections

from Bio             import SeqIO           as _SeqIO
from Bio.Seq         import Seq             as _Seq
from Bio.SeqRecord   import SeqRecord       as _SeqRecord
from Bio.SeqFeature  import SeqFeature      as _SeqFeature
from Bio.SeqFeature  import FeatureLocation as _FeatureLocation

from pygenome._pretty import pretty_str as _ps

from pydna.readers    import read as _read
from pydna.dseqrecord import Dseqrecord as _Dseqrecord
from pydna.amplify    import pcr        as _pcr
from pydna.assembly   import Assembly   as _Assembly
from pydna.primer     import Primer as _Primer

from pkg_resources   import resource_filename as _resource_filename

_pFA6a_kanMX4         = _read( _resource_filename("pygenome", "pFA6a-kanMX4.gb") )
_pFA6a_GFPS65T_kanMX6 = _read( _resource_filename("pygenome", "pFA6a-GFPS65T-kanMX6.gb") )

from pygenome.systematic_name import _systematic_name
from pygenome.standard_name   import _standard_name
from pygenome._data     import _data_files
from pygenome.intergenic     import  intergenic_sequence

data_dir = _os.path.join( _os.getenv("pygenome_data_dir"), "Saccharomyces_cerevisiae")

if _sys.version_info >= (3,8):
    import pickle  as _pickle
else:
    import pickle5 as _pickle

_feature_list = _pickle.load( open(_os.path.join(data_dir, "feature_list.pickle"), "rb" ) )
_not_deleted =_pickle.load(open(_os.path.join(data_dir, "not_done.pickle"), "rb" ))
_descriptions = _pickle.load( open(_os.path.join(data_dir, "systematic_to_description.pickle"), "rb" ) )
_primers = _pickle.load( open(_os.path.join(data_dir, "primers.pickle"), "rb" ) )

class Gene():

    def __init__(self, gene_name):
        self.sys = _ps(_systematic_name(gene_name))
        self.std = _ps(_standard_name(gene_name) or self.sys)
        self.chr = _ps(_data_files[ord(self.sys[1])-65])
        link = "<a href='http://www.yeastgenome.org/locus/{gene}' target='_blank'>{text}</a>"
        self.sgd_link = _ps(link.format(gene=self.sys, text="Gene {}/{}".format(self.std, self.sys)))

    def __repr__(self):
        return "Gene {}/{}".format(self.std, self.sys)

    def _repr_pretty_(self, p, cycle):
        '''returns a short string representation of the object'''
        p.text("Gene {}/{}".format(self.std, self.sys))

    def _repr_html_(self):
        return self.sgd_link

    def __len__(self):
        return len(self.locus(upstream=0, downstream=0))

    def efetch_link(self):
        with open(_os.path.join(data_dir, _data_files[ord(self.sys[1])-65])) as f:
            header = f.read(4*79)
        for line in header.splitlines():
            if line.startswith("VERSION"):
                key, accession = line.split()
                break
        _krom = _SeqIO.read(_os.path.join(data_dir, _data_files[ord(self.sys[1])-65]), "gb")
        feature = {f.qualifiers['locus_tag'][0] :  f for f in [f for f in _krom.features if f.type=="CDS"]}[self.sys]
        start = feature.location.start + 1
        stop  = feature.location.end
        strand = feature.location.strand
        return _ps(f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={accession}&strand={strand}&seq_start={start}&seq_stop={stop}&rettype=gb&retmode=text")

    def locus(self, upstream=1000, downstream=1000):
        _krom = _SeqIO.read(_os.path.join(data_dir, _data_files[ord(self.sys[1])-65]), "gb")
        feature = {f.qualifiers['locus_tag'][0] :  f for f in [f for f in _krom.features if f.type=="CDS"]}[self.sys]

        color = '#%02x%02x%02x' % (int(_random.uniform(150,255)),
                                   int(_random.uniform(150,255)),
                                   int(_random.uniform(150,255)),)

        feature.qualifiers.update({"ApEinfo_fwdcolor" : color,
                                   "ApEinfo_revcolor" : color,
                                   })

        start, stop = feature.location.start, feature.location.end
        #lcs = feature.extract(_krom)

        if self.sys[6]=="W":
            lcs = _krom[start-upstream:stop+downstream]
            lcs.description = "{} REGION: {}..{}".format(_krom.name,start-upstream+1, stop+downstream)
        else:
            lcs = _krom[start-upstream:stop+downstream].reverse_complement()
            lcs.description = "{} REGION: complement({}..{})".format(_krom.name, stop+downstream+1, start-upstream)

        lcs.name = _krom.name

        lcs.id = _krom.id

        #lcs.annotations

        lcs.pydna_code = _ps( "from pydna.genbank import Genbank\n"
                              "gb = Genbank('my@email.com')\n"
                              "seq = gb.nucleotide('{}')".format(lcs.description))
        return _Dseqrecord.from_SeqRecord(lcs,linear=True)


    def orf(self):
        '''Returns the open reading frame associated with a standard name
        (eg. CYC1) or a systematic name (eg. YJR048W).

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
        >>> from pygenome import saccharomyces_cerevisiae as sg
        >>> sg.stdgenes["TDH3"].cds()
        Dseqrecord(-999)
        >>> len(sg.stdgenes["TDH3"].cds())
        999
        '''

        return self.locus(upstream=0, downstream=0)


    def cds(self):
        '''

        '''
        ft = [f for f in self.orf().features if f.type=="CDS"]
        ft = ft.pop()
        return ft.extract(self.orf())


    def upstream_gene(self):
        '''Returns the coding sequence (cds) assciated with the gene upstream
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


    def downstream_gene(self):
        '''Returns the coding sequence (cds) assciated with the gene downstream
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


    def short_description(self):
        return _descriptions[self.sys]


    def promoter(self):
        '''Returns the sequence of the promoter assciated with
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
        >>> from pygenome import saccharomyces_cerevisiae as sg
        >>> sg.stdgenes["TDH3"].promoter()
        Dseqrecord(-698)
        >>> str(sg.stdgenes["TDH3"].promoter) == str(sg.sysgenes["YGR192C"].promoter)
        True
        >>>
        '''

        pr = intergenic_sequence(self.upstream_gene().sys, self.sys)
        pr.id = self.upstream_gene().sys+"_"+self.sys
        pr.name = "."
        #pr.description = "Intergenic sequence between upstream gene {} and downstream gene {}".format(self.upstream_gene().sys, self.sys)
        pr.features.append(_SeqFeature(_FeatureLocation(0, len(pr)),
                                      type = "promoter",
                                      strand = 1,
                                      qualifiers = {"note"              : "tp {} {}".format(self.upstream_gene().sys, self.sys),
                                                    "ApEinfo_fwdcolor": "#b1e6cc",
                                                    "ApEinfo_revcolor": "#b1e681",
                                                    }))
        return pr




    def terminator(self):
        '''Returns the sequence of the terminator assciated with
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
        terminator_pydna_code

        Examples
        --------
        >>> from pygenome import saccharomyces_cerevisiae as sg
        >>> sg.stdgenes["TDH3"].terminator()
        Dseqrecord(-580)
        >>> len(sg.stdgenes["TDH3"].terminator())
        580
        >>> sg.stdgenes["TDH3"].upstream_gene()
        Gene PDX1/YGR193C
        >>> str(sg.stdgenes["PDX1"].terminator().seq) == str(sg.stdgenes["TDH3"].promoter().seq)
        True
        '''


        tm = intergenic_sequence( self.sys, self.downstream_gene().sys )
        tm.id = self.sys + "_" + self.upstream_gene().sys
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


    def tandem(self):
        '''Returns True if a gene is expressed in the same direction on the chromosome as
        the gene immediatelly upstream.

        ::
                Tandem genes:

                Gene1 Gene2
                ----> ---->

                bidirectional genes:

                Gene1 Gene2
                ----> <----

        See also the bidirectional property.

        Returns
        -------
        out : bool
            Boolean; True or False

        '''
        return self.sys[6] == self.upstream_gene().sys[6]



    def bidirectional(self):
        '''Returns True if a gene is NOT expressed in the same direction on the chromosome as
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

        See also the tandem property.

        Returns
        -------
        out : bool
            Boolean; True or False
        '''
        return not self.tandem()

    def deletion_cassettes(self):

        if _not_deleted[self.sys]:
            warn("Deletion not done for this gene.")

        cassettes = []

        for p in _primers[self.sys]:

            upt = _Primer(p.UPTAG_primer_sequence)
            dnt = _Primer(p.DNTAG_primer_sequence)
            ups = _Primer(p.UPstream45_primer_sequence)
            dns = _Primer(p.DNstream45_primer_sequence)


            upt.id  = "UPTAG_primer_{}".format(self.sys)
            dnt.id  = "DNTAG_primer_{}".format(self.sys)

            ups.id  = "UPstream45_{}".format(self.sys)
            dns.id  = "DNstream45_{}".format(self.sys)

            if p.UPTAG_primer_sequence and p.DNTAG_primer_sequence:
                try:
                    inner_cassette = _pcr( upt, dnt, _pFA6a_kanMX4)
                except ValueError:
                    warn("PCR error in first cassette")
                    inner_cassette = None
                try:
                    outer_cassette = _pcr( ups, dns, inner_cassette)
                except ValueError:
                    warn("PCR error in second cassette")
                    outer_cassette = None
            else:
                try:
                    inner_cassette = _pcr( upt, dns, _pFA6a_kanMX4)
                except ValueError:
                    warn("PCR error in first cassette")
                    inner_cassette = None
                try:
                    outer_cassette = _pcr( ups, dns, inner_cassette)
                except ValueError:
                    warn("PCR error in second cassette")
                    outer_cassette = None

            if outer_cassette:
                outer_cassette.add_feature(0,len(outer_cassette), )

            cas = _collections.namedtuple('Cassette', ['outer_cassette',
                                                      'inner_cassette',
                                                      'UPTAG',
                                                      'DNTAG',
                                                      'UPstream45',
                                                      'DNstream45'])
            cassettes. append(cas(outer_cassette, inner_cassette, upt, dnt, ups, dns))

        return cassettes


    def deletion_loci(self):

        outer_cassettes = [x.outer_cassette for x in self.deletion_cassettes() if x.outer_cassette]

        loci = []

        for cas in outer_cassettes:

            locus = self.locus(upstream=1000, downstream=1000)

            asm = _Assembly( (locus, cas, locus))

            candidates = asm.assemble_linear()

            if candidates:
                kanmx4_gene = candidates[0]

                kanmx4_gene.name = "{}::KanMX4".format(self.sys.lower())

                kanmx4_gene.id = "{} locus with {} bp up and {} bp downstream DNA".format(kanmx4_gene.name, 1000, 1000)
            else:
                kanmx4_gene = None
                warn("No integration of cassette.")

            loci.append(kanmx4_gene)

        return loci


    def gfp_cassette(self):

        # We used the ‘Promoter’ program (courtesy of Joe DeRisi, publicly available at:
        # http://derisilab.ucsf. edu) to extract the last 40 nucleotides (excluding the
        # stop codon) of each ORF, as well as 40 nucleotides of genomic sequence
        # immediately following the stop codon of each ORF.
        #
        # We added the constant forward sequence from
        # the ‘Pringle’ oligonucleotide-directed homologous
        # recombination system (Longtine et al., 1998) to
        # the last 40 nucleotides of each ORF to create
        # the F2 oligo sequence, and the reverse complement
        # of the 40 nucleotides following each ORF
        # to the constant reverse sequence to create the R1
        # oligo sequence.
        #
        # last40wostop - CGGATCCCCGGGTTAATTAA F2
        # rc of 40 nucleotides following each - GAATTCGAGCTCGTTTAAAC R1
        #
        #
        # Longtine, M.S., McKenzie, A., 3rd, Demarini, D.J., Shah, N.G., Wach, A.,
        # Brachat, A., Philippsen, P., Pringle, J.R., 1998. Additional modules for
        # versatile and economical PCR-based gene deletion and modification in
        # Saccharomyces cerevisiae. Yeast 14, 953–961.
        #
        # AJ002682 pFA6a-GFPS65T-kanMX6

        F2 = "CGGATCCCCGGGTTAATTAA"
        R1 = "GAATTCGAGCTCGTTTAAAC"

        forward_tag_primer = self.cds()[-43:-3] + F2

        forward_tag_primer.id = "forward_tag_primer"

        reverse_tag_primer = self.terminator()[:40].reverse_complement() + R1

        reverse_tag_primer.id = "reverse_tag_primer"

        cassette = _pcr(forward_tag_primer, reverse_tag_primer, _pFA6a_GFPS65T_kanMX6)

        genome_gene_locus = _Dseqrecord( self.locus() )

        asm = _Assembly((genome_gene_locus, cassette, genome_gene_locus))

        candidate = asm.assemble_linear()[0]

        candidate.figure()

        cassette_in_genome = candidate

        GFP_fp = cassette_in_genome[1000:len(cassette_in_genome)-1000-(len(cassette_in_genome)-1000)%3+1].seq.translate( to_stop=True)

        assert (self.cds()[:-3]+_pFA6a_GFPS65T_kanMX6[41:776]).seq.translate(stop_symbol='') == GFP_fp

        return cassette
