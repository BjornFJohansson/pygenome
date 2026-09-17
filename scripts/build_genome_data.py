#!/usr/bin/env python3
"""Build sequence-free S288C metadata from sixteen local GenBank files."""
import argparse
from collections import defaultdict
import hashlib
import json
from pathlib import Path
import re

from Bio import SeqIO
from Bio.SeqFeature import ExactPosition


def location(feature, length):
    loc = feature.location
    if loc is None or loc.strand not in (-1, 1):
        raise ValueError(f"Unsupported strand/location: {loc}")
    if getattr(loc, "operator", "join") != "join":
        raise ValueError(f"Unsupported location operator: {loc}")
    intervals = []
    for part in loc.parts:
        if (type(part.start) is not ExactPosition or type(part.end) is not ExactPosition
                or part.ref or part.ref_db or part.strand != loc.strand):
            raise ValueError(f"Unsupported location: {loc}")
        start, end = int(part.start) + 1, int(part.end)
        if not 1 <= start <= end <= length:
            raise ValueError(f"Location outside chromosome bounds: {loc}")
        intervals.append([start, end])
    intervals.sort()
    if any(a[1] >= b[0] for a, b in zip(intervals, intervals[1:])):
        raise ValueError(f"Overlapping location parts: {loc}")
    return intervals


def build(paths):
    """Parse the supplied chromosomes; fail rather than drop invalid annotations."""
    genes, chromosomes, provenance = {}, {}, []
    aliases = defaultdict(set)
    cds_features, issues, conventions = [], [], []
    counts = {"gene": 0, "CDS": 0}
    for number, path in enumerate(sorted(map(Path, paths)), 1):
        record = SeqIO.read(path, "genbank")
        accession, length = record.id, len(record)
        if not re.fullmatch(r"[A-Za-z0-9_]+\.\d+", accession):
            issues.append(f"{path.name}: missing versioned accession: {accession}")
        if accession in chromosomes:
            issues.append(f"{path.name}: duplicate chromosome {accession}")
        chromosomes[accession] = {"length": length, "number": number}
        provenance.append({"filename": path.name, "accession": accession,
                           "sha256": hashlib.sha256(path.read_bytes()).hexdigest()})
        for feature in record.features:
            if feature.type not in counts:
                continue
            counts[feature.type] += 1
            tags = feature.qualifiers.get("locus_tag", [])
            if len(tags) != 1 or not tags[0]:
                issues.append(f"{path.name}: {feature.type} missing/ambiguous locus_tag")
                continue
            tag = tags[0]
            try:
                intervals = location(feature, length)
            except ValueError as exc:
                issues.append(f"{path.name}: {tag} {feature.type}: {exc}")
                continue
            if feature.type == "CDS":
                cds_features.append((tag, accession, feature.location.strand, intervals))
                continue
            if tag in genes:
                issues.append(f"{path.name}: duplicate/conflicting gene {tag}")
                continue
            genes[tag] = {"accession": accession, "strand": feature.location.strand,
                          "gene": intervals, "cds": []}
            for name in [tag, *feature.qualifiers.get("gene", []),
                         *feature.qualifiers.get("gene_synonym", [])]:
                aliases[name.strip().upper()].add(tag)
            match = re.fullmatch(r"Y([A-P])[LR]\d{3}([WC])(?:-[A-Z])?", tag)
            if match and (ord(match[1]) - ord("A") + 1 != number
                          or (1 if match[2] == "W" else -1) != feature.location.strand):
                conventions.append(tag)
    for tag, accession, strand, intervals in cds_features:
        if tag not in genes:
            issues.append(f"Unmatched CDS: {tag}")
            continue
        gene = genes[tag]
        if (gene["accession"] != accession or gene["strand"] != strand
                or any(not gene["gene"][0][0] <= a <= b <= gene["gene"][-1][1]
                       for a, b in intervals)):
            issues.append(f"CDS inconsistent with gene: {tag}")
        elif intervals not in gene["cds"]:
            gene["cds"].append(intervals)
    if issues:
        raise ValueError("Invalid source annotations:\n" + "\n".join(issues))
    for gene in genes.values():
        gene["cds"].sort()
    return {"schema_version": 1, "genome": "S288C", "chromosomes": chromosomes,
            "genes": genes, "aliases": {k: sorted(v) for k, v in aliases.items()},
            "provenance": provenance,
            "audit": {"feature_counts": counts, "name_convention_mismatches": sorted(conventions),
                      "multiple_cds": sorted(k for k, v in genes.items() if len(v["cds"]) > 1),
                      "genes_without_cds": sum(not v["cds"] for v in genes.values())}}


def serialize(data):
    return json.dumps(data, indent=2, sort_keys=True) + "\n"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    paths = [args.input_dir / f"chr{i:02}.gbf" for i in range(1, 17)]
    missing = [str(p) for p in paths if not p.is_file()]
    if missing:
        parser.error("Missing chromosome files: " + ", ".join(missing))
    try:
        data = build(paths)
    except ValueError as exc:
        parser.exit(1, str(exc) + "\n")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(serialize(data), encoding="utf-8")
    print(json.dumps(data["audit"], indent=2))
    collisions = {k: v for k, v in data["aliases"].items() if len(v) > 1}
    print(f"Retained {len(collisions)} ambiguous aliases; see generated aliases mapping.")


if __name__ == "__main__":
    main()
