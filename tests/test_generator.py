from pathlib import Path
import json
import pytest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, SimpleLocation, CompoundLocation, BeforePosition
from scripts.build_genome_data import build, serialize

ROOT = Path(__file__).resolve().parents[1]


def feature(kind, start, end, tag="A", strand=1, **qualifiers):
    return SeqFeature(SimpleLocation(start, end, strand=strand), type=kind,
                      qualifiers={"locus_tag": [tag], **qualifiers})


def write(tmp_path, features):
    record = SeqRecord(Seq("A" * 100), id="TEST.1", name="TEST", annotations={"molecule_type": "DNA"})
    record.features = features
    path = tmp_path / "chr01.gbf"
    SeqIO.write(record, path, "genbank")
    return path


def test_fixture(tmp_path):
    gene = feature("gene", 9, 40, strand=-1, gene=["Named"], gene_synonym=["one", "two"])
    cds = feature("CDS", 9, 20, strand=-1)
    cds.location = CompoundLocation([SimpleLocation(29, 40, -1), SimpleLocation(9, 20, -1)])
    path = write(tmp_path, [cds, gene, feature("gene", 59, 70, "B", gene_synonym=["one"])])
    data = build([path])
    assert data["genes"]["A"]["gene"] == [[10, 40]]
    assert data["genes"]["A"]["cds"] == [[[10, 20], [30, 40]]]
    assert data["aliases"]["ONE"] == ["A", "B"]
    assert data["aliases"]["TWO"] == ["A"]
    assert data["genes"]["B"]["cds"] == []
    assert serialize(data) == serialize(build([path]))


@pytest.mark.parametrize("case", ["missing", "unmatched", "duplicate", "fuzzy", "order", "bounds", "strand"])
def test_invalid_annotations(tmp_path, case):
    genes = [feature("gene", 10, 30)]
    if case == "missing":
        genes[0].qualifiers = {}
    elif case == "unmatched":
        genes += [feature("CDS", 10, 30, "B")]
    elif case == "duplicate":
        genes += [feature("gene", 40, 50)]
    elif case == "fuzzy":
        genes[0].location = SimpleLocation(BeforePosition(10), 30, 1)
    elif case == "order":
        genes[0].location = CompoundLocation([SimpleLocation(10, 20, 1), SimpleLocation(25, 30, 1)], "order")
    elif case == "bounds":
        genes[0].location = SimpleLocation(10, 110, 1)
    elif case == "strand":
        genes += [feature("CDS", 10, 30, strand=-1)]
    with pytest.raises(ValueError):
        build([write(tmp_path, genes)])


def test_multiple_cds_and_deduplication(tmp_path):
    data = build([write(tmp_path, [feature("gene", 10, 30), feature("CDS", 10, 30),
                                  feature("CDS", 10, 30), feature("CDS", 15, 30)])])
    assert data["genes"]["A"]["cds"] == [[[11, 30]], [[16, 30]]]
    assert data["audit"]["multiple_cds"] == ["A"]


def test_complete_dataset():
    paths = [ROOT / "data" / f"chr{i:02}.gbf" for i in range(1, 17)]
    if not all(p.exists() for p in paths):
        pytest.skip("Source GenBank files are intentionally excluded from distributions")
    data = build(paths)
    assert len(data["genes"]) == data["audit"]["feature_counts"]["gene"] == 6424
    assert data["audit"]["feature_counts"]["CDS"] == 6001
    assert data["audit"]["genes_without_cds"] == 423
    assert not data["audit"]["name_convention_mismatches"]
    assert serialize(data) == (ROOT / "src/pygenome/S288C.json").read_text()
    for entry in data["genes"].values():
        length = data["chromosomes"][entry["accession"]]["length"]
        for intervals in [entry["gene"], *entry["cds"]]:
            assert all(1 <= a <= b <= length for a, b in intervals)
