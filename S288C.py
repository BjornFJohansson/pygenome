#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Access to the Saccharomyces cerevisiae genome from Python.

Sequences can be accessed as Biopython Bio.SeqRecord or pydna.Dseqrecord
objects.

http://sgd-archive.yeastgenome.org/sequence/S288C_reference/NCBI_genome_source

"""
import os as _os
from email.utils import parsedate_to_datetime as _ptod
from tqdm import tqdm as _tqdm
import requests as _requests
from pathlib import Path as _Path
import sys as _sys
from urllib.parse import urlparse as _up
from os.path import basename as _bn
from Bio import SeqIO as _SeqIO
from pydna.genbankrecord import GenbankRecord as _gbr
from pyfaidx import Fasta as _Fasta
from send2trash import send2trash as _send2trash

if _sys.version_info >= (3, 8):
    import pickle as _pickle
else:
    import pickle5 as _pickle  # pragma: no cover
from pydna.primer import Primer as _Primer
from pydna.amplify import pcr as _pcr
from pydna.assembly import Assembly as _Assembly
from pydna.readers import read as _read
from warnings import warn as _warn
from collections import namedtuple as _ntup
import appdirs as _appdirs
import pickletools as _pickletools
from pkg_resources import resource_filename as _resource_filename
from zipfile import ZipFile as _zf
from calendar import timegm as _timegm
from datetime import timedelta as _timedelta
pygenome_data_dir = _Path(_appdirs.user_data_dir("pygenome"))

# /home/bjorn/.local/share/pygenome/S288C
data_dir = pygenome_data_dir/"S288C"

# set the chromosome_urls string variable
try:
    code = (data_dir/"chromosome_urls.py").read_text()
except FileNotFoundError:
    chromosome_urls = ""
else:
    exec(code, globals(), locals())
    assert chromosome_urls


class Gene():
    """docstring."""

    def __init__(self,
                 sysname,
                 stdname,
                 start,
                 end,
                 strand,
                 accession,
                 fa_pth,
                 fe_pth):
        self.sysname = sysname
        self.stdname = stdname
        # assert start <= end
        self.start = start
        self.end = end
        self.strand = strand
        self.accession = accession
        self.fa_pth = fa_pth
        self.fe_pth = fe_pth
        self.pred = None
        self.succ = None

    def __repr__(self):
        """docstring."""
        return "Gene {}/{}".format(self.stdname, self.sysname)

    def _repr_pretty_(self, p, cycle):
        """Short string representation of the object."""
        p.text(self.__repr__())

    def _repr_html_(self):
        """docstring."""
        link = (f"<a href='http://www.yeastgenome.org/locus/{self.sysname}' "
                f"target='_blank'>{self.__repr__()}</a>")
        return link

    def __len__(self):
        """docstring."""
        return self.end - self.start

    def _slice(self, start, end, strand):
        s = _Fasta(str(self.fa_pth))[self.accession][start:end]
        result = _gbr.from_string(str(s if strand == 1 else -s),
                                  item=self.accession,
                                  start=start+1,
                                  stop=end,
                                  strand={1: 1, -1: 2}[strand])
        features = _pickle.load(self.fe_pth.open("rb"))
        if strand == 1:
            result.features = [f._shift(-start) for f in features if
                               f.location.start >= start
                               and f.location.end <= end]
            result.description = f"{self.accession} REGION: {start+1}..{end}"
        else:
            result.features = [f._flip(end) for f in features if
                               f.location.start >= start
                               and f.location.end <= end]
            result.description = (f"{self.accession} REGION: "
                                  f"complement({start+1}..{end})")
        return result

    @property
    def dna(self):
        """docstring."""
        return self._slice(self.start, self.end, self.strand)

    @property
    def cds(self):
        """docstring."""
        # There should be one and only one CDS feature
        ft, = [ft for ft in self.dna.features if ft.type=="CDS"]
        return ft.extract(self.dna)

    @property
    def promoter(self):
        """docstring."""
        if self.strand == 1:
            result = self._slice(self.pred.end, self.start, 1)
        else:
            result = self._slice(self.end, self.succ.start, -1)
        return result

    @property
    def terminator(self):
        """docstring."""
        if self.strand == 1:
            result = self._slice(self.end, self.succ.start, 1)
        else:
            result = self._slice(self.pred.end, self.start, -1)
        return result

    @property
    def transcriptional_unit(self):
        """docstring."""
        return self._slice(self.pred.end, self.succ.start, self.strand)

    tu = transcriptional_unit

    def locus(self, upstream=1000, downstream=1000):
        """docstring."""
        locus = self._slice(max(self.start-upstream, 0),
                           self.end+downstream,
                           self.strand)


        return locus

    @property
    def tandem(self):
        """docstring."""
        return {1:  self.pred.strand == 1,
                -1: self.succ.strand == -1}[self.strand]

    @property
    def divergent(self):
        """docstring."""
        return not(self.tandem)

    @property
    def deletion_cassettes(self):
        """docstring."""
        primers = _pickle.load(
            (data_dir/"Deletion_primers_PCR_sizes.pickle").open("rb"))

        cassettes = []

        pFA6a_kanMX4 = _pickle.load(
            (data_dir/"pFA6a-kanMX4.pickle").open("rb"))

        for p in primers[self.sysname]:

            upt = _Primer(p[8])
            dnt = _Primer(p[9])
            ups45 = _Primer(p[10])
            dns45 = _Primer(p[11])

            upt.id = f"UPTAG_{self.sysname}"
            dnt.id = f"DNTAG_{self.sysname}"
            ups45.id = f"UPstream45_{self.sysname}"
            dns45.id = f"DNstream45_{self.sysname}"

            #    ch
            # A   1
            # E   5
            # M  13
            if upt and dnt:
                try:
                    inner_cassette = _pcr(upt, dnt, pFA6a_kanMX4)
                except ValueError:
                    _warn("PCR error in first cassette %s"
                          "\nprimer %s"
                          "\nprimer %s"
                          % (self.sysname,
                             upt.id,
                             dnt.id))
                    inner_cassette = None
                finally:
                    inner_cassette.add_feature(0, len(inner_cassette),)
                try:
                    outer_cassette = _pcr(ups45, dns45, inner_cassette)
                except ValueError:
                    _warn("PCR error in 2nd cassette %s"
                          "\nprimer %s"
                          "\nprimer %s"
                          % (self.sysname,
                             upt.id,
                             dnt.id))
                    outer_cassette = None
                finally:
                   outer_cassette.add_feature(0, len(outer_cassette),)
            else:
                try:
                    inner_cassette = _pcr(upt, dns45, pFA6a_kanMX4)
                except ValueError:
                    _warn("PCR error in first cassette %s"
                          "\nprimer %s"
                          "\nprimer %s"
                          % (self.sysname,
                             upt.id,
                             dns45.id))
                    inner_cassette = None
                finally:
                    inner_cassette.add_feature(0, len(inner_cassette),)
                try:
                    outer_cassette = _pcr(ups45, dns45, inner_cassette)
                except ValueError:
                    _warn("PCR error in 2nd cassette %s"
                          "\nprimer %s"
                          "\nprimer %s"
                          % (self.sysname,
                             upt.id,
                             dns45.id))
                    outer_cassette = None
                finally:
                    outer_cassette.add_feature(0, len(outer_cassette),)

            cassettes.append(outer_cassette)

        return cassettes

    def cassette_integration_locus(self, cassette=None, locus=None, limit=25):
        """docstring."""
        if not cassette:
            cassette = self.deletion_cassettes[0]

        if not locus:
            locus = self.locus()

        asm = _Assembly((locus, cassette, locus), limit=limit)

        candidates = asm.assemble_linear()

        if candidates:

            candidate = candidates[0]

            candidate.name = f"{self.sysname.lower()}::{cassette.name}"

            candidate.id = (f"{candidate.name} locus with "
                            "1000 bp up and downstream DNA")
        else:
            candidate = None
            _warn("No integration of cassette.")

        return candidate


    @property
    def gfp_cassette_hismx(self):
        """
        GFP tagging cassette with HIS3MX6 marker.

        According to Huh 2003:

        Huh, Won-Ki, James V. Falvo, Luke C. Gerke, Adam S. Carroll,
        Russell W. Howson, Jonathan S. Weissman, and Erin K. O’Shea. 2003.
        “Global Analysis of Protein Localization in Budding Yeast.”
        *Nature* 425 (6959) (October): 686–691.

        https://yeastgfp.yeastgenome.org
        """

        tmpl = _pickle.load((data_dir/
                             "pFA6a-GFPS65T-HIS3MX6.pickle").open("rb"))

        GFPprimers = _pickle.load((data_dir/
                                   "yeastGFPOligoSequence.pickle").open("rb"))

        try:
            sys, std, fp, rp, cp = GFPprimers.get(self.sysname)
        except TypeError:
            _warn("No cassette for this gene.")
            return None
        else:
            fp = _Primer(fp)
            rp = _Primer(rp)
            cassette = _pcr(fp, rp, tmpl)

        cassette.name = f"{self.sysname}-GFP"

        return cassette



    @property
    def gfp_cassette_kanmx(self):
        """
        GFP tagging cassette with kanMX6 marker.

        According to Howson 2005:

        Howson, Russell, Won-Ki Huh, Sina Ghaemmaghami, James V. Falvo, Kiowa
        Bower, Archana Belle, Noah Dephoure, Dennis D. Wykoff,
        Jonathan S. Weissman, and Erin K. O’Shea. 2005.
        “Construction, Verification and Experimental Use of Two Epitope-Tagged
        Collections of Budding Yeast Strains.” *Comparative and Functional
        Genomics* 6 (1-2): 2–16.

        Briefly from the materials section:

        We used the ‘Promoter’ program (courtesy of Joe DeRisi, publicly
        available at: http://derisilab.ucsf.edu) to extract the last
        40 nucleotides (excluding the stop codon) of each ORF, as well as
        40 nucleotides of genomic sequence immediately following the stop
        codon of each ORF.

        We added the constant forward sequence from the ‘Pringle’
        oligonucleotide-directed homologous recombination system
        (Longtine et al., 1998) to the last 40 nucleotides of each ORF to
        create the F2 oligo sequence, and the reverse complement of the 40
        nucleotides following each ORF to the constant reverse sequence to
        create the R1 oligo sequence.

        >F2
        40 nt upstream stop codon - CGGATCCCCGGGTTAATTAA

        >R1
        rc of 40 nt downstream - GAATTCGAGCTCGTTTAAAC

        Template was lasmid pFA6a-GFPS65T-kanMX6 (AJ002682.1)
        """
        F2 = "CGGATCCCCGGGTTAATTAA"
        R1 = "GAATTCGAGCTCGTTTAAAC"

        tmpl = _pickle.load((data_dir/"pFA6a-GFPS65T-kanMX6.pickle").open("rb"))

        from pydna.primer import Primer

        forward_tag_primer = Primer(self.cds[-43:-3] + F2)

        forward_tag_primer.id = "forward_tag_primer"

        reverse_tag_primer = Primer(str(self.terminator.seq[:40].reverse_complement()) + R1)

        reverse_tag_primer.id = "reverse_tag_primer"

        cassette = _pcr(forward_tag_primer,
                        reverse_tag_primer,
                        tmpl)

        cassette.name = f"{self.sysname}-GFP"

        return cassette


    @property
    def short_description(self):
        """docstring."""
        features = _pickle.load(self.fe_pth.open("rb"))
        features = [f for f in features if
                    f.location.start >= self.start
                    and f.location.end <= self.end]
        return (features[-1].qualifiers.get("note") or [""])[0]


def extract_data():
    """docstring."""
    zf = _zf(_resource_filename("pygenome",
                                "saccharomyces_cerevisiae/S288C.zip"), "r")
    zf.extractall(data_dir)
    for member in zf.infolist():
        date_time = _timegm(member.date_time)
        _os.utime(data_dir/member.filename,
                  (date_time, date_time))


def check_data_files():
    """docstring."""
    zf = _zf(_resource_filename("pygenome",
                                "saccharomyces_cerevisiae/S288C.zip"), "r")
    for url in chromosome_urls.splitlines():
        response = _requests.get(url, stream=True)
        remote_time_stamp = _ptod(response.headers.get('last-modified') or 0
                                  ).timestamp()
        fa_pth = _Path(_bn(_up(url).path)).with_suffix(".fasta")
        local_time_stamp = _timegm(zf.getinfo("chr01.fasta").date_time)

        if local_time_stamp < remote_time_stamp:
            print(f"{fa_pth.stem} is out of date.")
        else:
            print(f"{fa_pth.stem} is up to date")
    zf.close()


def update_data_files():
    """docstring."""
    for url in chromosome_urls.splitlines():

        response = _requests.get(url, stream=True)
        remote_time_stamp = (_ptod((
            response.headers.get('last-modified') or 0)) + _timedelta(seconds=5)).timestamp()
        total = int(response.headers.get('content-length'))

        fn = _bn(_up(url).path)
        gb_pth = data_dir/fn

        with open(gb_pth, 'wb') as f:
            for data in _tqdm(response.iter_content(),
                              desc=gb_pth.name,
                              total=total):
                f.write(data)
        krom = _SeqIO.read(gb_pth, "gb")

        fa_pth = gb_pth.with_suffix(".fasta")

        _SeqIO.write(krom, fa_pth, "fasta")
        _os.utime(fa_pth, times=(remote_time_stamp,)*2)
        pkl = _pickle.dumps(tuple(krom.features), 5)
        pkl = _pickletools.optimize(pkl)

        fe_pth = fa_pth.with_name(f"{fa_pth.stem}_feature_tuple.pickle")

        with fe_pth.open("wb") as f:
            f.write(pkl)
        _os.utime(fe_pth, times=(remote_time_stamp,)*2)

        try:
            _send2trash(gb_pth)
        except OSError:
            pass






def gene_dicts():
    """docstring."""
    genedict = {}
    for url in chromosome_urls.splitlines():
        fn = _bn(_up(url).path)
        gb_pth = data_dir/fn
        fa_pth = gb_pth.with_suffix(".fasta")
        with fa_pth.open("r") as f:
            accession = f.readline().lstrip(">").split().pop(0)
        fe_pth = fa_pth.with_name(f"{fa_pth.stem}_feature_tuple.pickle")
        genelist = []
        features = _pickle.load(fe_pth.open("rb"))
        for f in features:
            if f.type == 'CDS':
                sysname = f.qualifiers["locus_tag"][0]
                stdname = (f.qualifiers.get("gene") or [None]).pop()
                genelist.append(Gene(sysname,
                                     stdname,
                                     int(f.location.start),
                                     int(f.location.end),
                                     f.location.strand,
                                     accession,
                                     fa_pth,
                                     fe_pth))
        for i in range(1, len(genelist)-1):
            genelist[i].pred = genelist[i-1]
            genelist[i].succ = genelist[i+1]

        genelist[0].succ = genelist[1]
        genelist[-1].pred = genelist[-2]

        genedict.update({g.sysname: g for g in genelist})

    stdgenedict = {}

    for key, val in genedict.items():
        stdname = val.stdname
        if stdname:
            stdgenedict.update(((stdname, val),))

    return genedict, stdgenedict


def _feature_list(featuretype=""):
    """docstring."""
    flist = []
    for url in chromosome_urls.splitlines():
        fn = _bn(_up(url).path)
        fa_pth = (data_dir/fn).with_suffix(".fasta")
        ch = _read(fa_pth)
        fe_pth = fa_pth.with_name(f"{fa_pth.stem}_feature_tuple.pickle")
        ch.features = _pickle.load(fe_pth.open("rb"))
        for f in ch.features:
            if f.type != featuretype:
                continue
            flist.append(ch[f.location.start:f.location.end])
    return flist


def centromere_list():
    """docstring."""
    return _feature_list("centromere")


def rep_origin_list():
    """docstring."""
    return _feature_list("rep_origin")


def tRNA_list():
    """docstring."""
    return _feature_list("tRNA")


def genes_not_deleted_tuple():
    """docstring."""
    natuple = _pickle.load((data_dir/"ORFs_not_available.pickle").open("rb"))
    return natuple


if __name__ == "__main__":

    # from pygenome.saccharomyces_cerevisiae.S288C import gene_dicts

    sysgenes, stdgenes = gene_dicts()

    c = stdgenes["CYC1"]

    c.gfp_cassette

    cl = centromere_list()


    rl = rep_origin_list()


    tl = tRNA_list()


    """

    chr = g.promoter + g.cds + g.terminator
    chr.features
    [SeqFeature(FeatureLocation(ExactPosition(1029), ExactPosition(7185), strand=1), type='gene'),
     SeqFeature(FeatureLocation(ExactPosition(1029), ExactPosition(7185), strand=1), type='mRNA'),
     SeqFeature(FeatureLocation(ExactPosition(1029), ExactPosition(7185), strand=1), type='CDS')]
    chr.features
    [SeqFeature(FeatureLocation(ExactPosition(0), ExactPosition(1029), strand=1), type='promoter'),
     SeqFeature(FeatureLocation(ExactPosition(1029), ExactPosition(7185), strand=1), type='gene'),
     SeqFeature(FeatureLocation(ExactPosition(1029), ExactPosition(7185), strand=1), type='mRNA'),
     SeqFeature(FeatureLocation(ExactPosition(1029), ExactPosition(7185), strand=1), type='CDS'),
     SeqFeature(FeatureLocation(ExactPosition(7185), ExactPosition(7674), strand=1), type='terminator')]

    """

    # from pygenome.saccharomyces_cerevisiae.S288C import gene_dicts

    sysgenes, stdgenes = gene_dicts()

    # ProjectPygenome

    # <UTR1 ISY1>    ISY1 and UTR1 are both divergent
    # ISY1> OSM1>    OSM1 is tandem
    # <RAD23 <ANP1   RAD23 is tandem

    i = stdgenes["ISY1"]
    u = stdgenes["UTR1"]
    o = stdgenes["OSM1"]
    s = stdgenes["SNC1"]
    r = stdgenes["RAD23"]
    a = stdgenes["ANP1"]

    assert o.tandem == True
    assert r.tandem == True
    assert i.tandem == False
    assert u.tandem == False
    assert i.divergent == True
    assert u.divergent == True

    assert i.promoter.seq == u.promoter.seq.rc()
    assert i.terminator.seq == o.promoter.seq
    assert a.terminator.seq == r.promoter.seq

    for g in [i, u, o, s, r, a]:
        assert (g.promoter + g.dna + g.terminator).seq == g.tu.seq

    y = stdgenes["YRA1"]
    t = stdgenes["TAD3"]





    c = stdgenes["CYC1"]

    c.dna

    from pydna.genbank import Genbank
    gb = Genbank('bjornjobb@gmail.com')
    seq = gb.nucleotide('BK006943.2',
                        seq_start=526335,
                        seq_stop=526664,
                        strand=1)

    from pydna.genbank import genbank

    cc = genbank("BK006943.2 526335-526664")

    assert c.dna.seq == cc.seq == seq.seq

    prom = c.promoter

    item = repr(c.promoter).strip(")").split(" ", maxsplit=1)[1]

    prom_from_repr = genbank(item)

    prom.pydna_code()

    from pydna.genbank import Genbank
    gb = Genbank('bjornjobb@gmail.com')
    seq = gb.nucleotide('BK006943.2',
                        seq_start=525382,
                        seq_stop=526334,
                        strand=1)

    assert prom.seq==prom_from_repr.seq==seq.seq

    term = c.terminator

    item = repr(c.terminator).strip(")").split(" ", maxsplit=1)[1]

    term_from_repr = genbank(item)

    term.pydna_code()
    from pydna.genbank import Genbank
    gb = Genbank('bjornjobb@gmail.com')
    seq = gb.nucleotide('BK006943.2',
                        seq_start=526665,
                        seq_stop=526883,
                        strand=1)

    assert term.seq == term_from_repr.seq == seq.seq











"""

accession				    "BK006943.2"
dna
cds
locus
promoter
terminator
divergent
end                         526664
fa_pth                      PosixPath("/home/bjorn/.local/share/pygenome/S288C/chr10.fasta")
fe_pth                      PosixPath("/home/bjorn/.local/share/pygenome/S288C/chr10_feature_tuple.pickle")
pred                        Gene ANB1/YJR047C
short_description
start                       526334
stdname                     "CYC1"
strand                      1
succ                        Gene UTR1/YJR049C
sysname                     "YJR048W
tandem
transcriptional_unit / tu
deletion_cassettes
deletion_loci
gfp_cassette
gfp_tag_locus


accession				    "BK006945.2"
dna
cds
locus
promoter
terminator
divergent
end                         765265
fa_pth                      PosixPath('/home/bjorn/.local/share/pygenome/S288C/chr12.fasta')
fe_pth                      PosixPath('/home/bjorn/.local/share/pygenome/S288C/chr12_feature_tuple.pickle')
pred                        Gene NKP2/YLR315W
short_description
start                       766358
stdname                     "TAD3"
strand                      -1
succ                        Gene EST2/YLR318W
sysname                     "YLR316C"
tandem
transcriptional_unit / tu
deletion_cassettes
deletion_loci
gfp_cassette
gfp_tag_locus




stdname                     "TAD3"
sysname                     "YLR316C"
accession				    "BK006945.2"
strand                      -1
start                       766358
end                         765265
tandem                      False
divergent                   True
short_description           'Subunit of tRNA-specific adenosine-34 deaminase; forms a heterodimer with...'

dna                         Gbnk(-1093 BK006945.2 765266-766358)
cds                         Dseqrecord(-969)
locus                       Gbnk(-3093 BK006945.2 764266-767358)
promoter                    Gbnk(-183 BK006945.2 766359-766541)
terminator
transcriptional_unit / tu   Gbnk(-1272 BK006945.2 765270-766541)

fa_pth                      Path('/home/user/.local/share/pygenome/S288C/chr12.fasta')
fe_pth                      Path('/home/user/.local/share/pygenome/S288C/chr12_feature_tuple.pickle')

pred                        Gene NKP2/YLR315W
succ                        Gene EST2/YLR318W

deletion_cassettes
deletion_loci
gfp_cassette
gfp_tag_locus

"""
