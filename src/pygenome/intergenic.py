#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os            as _os
import itertools     as _itertools
from pydna.dseqrecord import Dseqrecord as _Dseqrecord

from Bio             import SeqIO as _SeqIO

from pygenome._pretty import pretty_str as _ps

from pygenome.systematic_name import _systematic_name

from pygenome._data import _data_files

data_dir = _os.path.join( _os.getenv("pygenome_data_dir"), "Saccharomyces_cerevisiae")

def intergenic_sequence(upgene, dngene):
    '''This function will return the intergenic sequence between two genes on the same chromosome.

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
    >>> from pygenome import saccharomyces_cerevisiae as sg
    >>> sg.stdgenes["TDH3"]
    Gene TDH3/YGR192C
    >>> sg.stdgenes["TDH3"].upstream_gene()
    Gene PDX1/YGR193C
    >>> sg.sysgenes["YGR193C"].upstream_gene()
    Gene XKS1/YGR194C
    >>> from pygenome.intergenic import intergenic_sequence
    >>> len(intergenic_sequence("YGR192C", "YGR193C"))
    698
    >>> len(intergenic_sequence("YGR192C", "YGR194C"))
    2262
    >>>

    '''

    upgene = _systematic_name(upgene)
    dngene = _systematic_name(dngene)

    if upgene[1] != dngene[1]:
        raise Exception("Both genes has to be on the same chromosome.")

    _krom = _SeqIO.read(_os.path.join(data_dir, _data_files[ord(upgene[1])-65]),"gb")
    cds  = {f.qualifiers['locus_tag'][0] :  f for f in [f for f in _krom.features if f.type=="CDS"]}
    upfeature = cds[upgene]
    startup, stopup  = upfeature.location.start,upfeature.location.end
    dnfeature = cds[dngene]
    startdn, stopdn = dnfeature.location.start, dnfeature.location.end

    assert sorted( (startup, stopup, startdn, stopdn) ) == list(_itertools.chain.from_iterable(sorted( [sorted((startup,stopup)),sorted((startdn, stopdn))] )))

    length,a,b = min([ (abs(a-b), a, b) for a, b in _itertools.product((startup,stopup),(startdn,stopdn))])

    if a<b:
        igs = _Dseqrecord.from_SeqRecord(_krom[a:b],linear=True)
        igs.description = "{} REGION: {}..{}".format(_krom.id, a+1,b)
    else:
        igs = _Dseqrecord.from_SeqRecord(_krom[b:a].reverse_complement(),linear=True)
        igs.description = "{} REGION: complement({}..{})".format(_krom.id, b+1,a)

    igs.name = _krom.name

    igs.id = _krom.id

    igs.pydna_code = _ps("gb = pydna.Genbank('my@email.com')\n"
                         "seq = gb.nucleotide('{}')".format(igs.description))
    return igs