"""
Extract a feature's sequence from a GenBank file using gb-io.

Targets a gene by locus_tag, e.g.:
    gene            complement(53260..54696)
                    /gene="ACT1"
                    /locus_tag="YFL039C"
"""

import gb_io

_COMPLEMENT_TABLE = bytes.maketrans(b"ACGTNacgtn", b"TGCANtgcan")

def reverse_complement(seq: bytes) -> bytes:
    return seq.translate(_COMPLEMENT_TABLE)[::-1]

records = gb_io.load("chr06.gbf") # list of gb_io.Record

for record in records:
    for feature in record.features:
        if feature.kind != "gene":
            continue
        qd = {}
        for q in feature.qualifiers:
            qd.setdefault(q.key, []).append(q.value)
        if qd.get("locus_tag") == ["YFL039C"]:
            break

start, stop = sorted((feature.location.start, feature.location.end))
sub = bytes(record.sequence[start - 1000: stop + 1000])

if isinstance(feature.location, gb_io.Complement):
    sub = reverse_complement(sub)
elif isinstance(feature.location, gb_io.Join):
    b"".join(bytes(record.sequence[part.location.start - 1 : part.location.end]) for part in feature.locations)

print(sub.decode())
