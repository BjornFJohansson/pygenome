import pytest
import pygenome as pg


def test_act1():
    assert pg.gene(" act1 ") == "BK006940.2 REGION: complement(53260..54696)"
    assert pg.gene("ABY1") == pg.gene("END7") == pg.gene("YFL039C") == pg.gene("ACT1")
    assert pg.cds("ACT1") == "BK006940.2 REGION: complement(join(53260..54377,54687..54696))"
    assert pg.locus("ACT1") == "BK006940.2 REGION: complement(52260..55696)"
    assert pg.locus("ACT1", upstream=10, downstream=20) == "BK006940.2 REGION: complement(53240..54706)"
    assert pg.promoter("ACT1") == pg.terminator("YPT1") == "BK006940.2 REGION: complement(54697..55365)"
    assert pg.cds("YPT1") == pg.gene("YPT1")


@pytest.fixture
def synthetic(monkeypatch):
    def entry(a, b, strand=1):
        return {"accession": "TEST.1", "strand": strand, "gene": [[a, b]], "cds": [[[a, b]]]}
    data = {"chromosomes": {"TEST.1": {"length": 100}},
            "genes": {"A": entry(10, 20), "B": entry(40, 50, -1), "C": entry(70, 80)},
            "aliases": {"A": ["A"], "B": ["B"], "C": ["C"], "AMB": ["A", "C"]}}
    monkeypatch.setattr(pg, "_data", lambda: data)
    return data


def test_gal_shared_promoter():
    assert pg.promoter("gal1") == "BK006936.2 REGION: 278353..279020"
    assert pg.promoter("gal10") == "BK006936.2 REGION: complement(278353..279020)"
    # The overlapping ncRNA remains available through ordinary gene lookup.
    assert pg.gene("YNCB0008W") == "BK006936.2 REGION: 276805..280645"


@pytest.mark.parametrize("interval", [[[1, 100]], [[25, 35]], [[55, 65]]])
def test_noncoding_neighbors_do_not_bound_gaps(synthetic, interval):
    synthetic["genes"]["RNA"] = {
        "accession": "TEST.1", "strand": 1, "gene": interval, "cds": []}
    assert pg.terminator("A") == "TEST.1 REGION: 21..39"
    assert pg.terminator("B") == "TEST.1 REGION: complement(21..39)"
    assert pg.promoter("B") == "TEST.1 REGION: complement(51..69)"
    assert pg.promoter("C") == "TEST.1 REGION: 51..69"


def test_synthetic_intervals(synthetic):
    assert pg.gene("A") == pg.cds("A") == "TEST.1 REGION: 10..20"
    assert pg.locus("A") == "TEST.1 REGION: 1..100"
    assert pg.locus("B") == "TEST.1 REGION: complement(1..100)"
    assert pg.locus("A", upstream=3, downstream=7) == "TEST.1 REGION: 7..27"
    assert pg.locus("B", upstream=3, downstream=7) == "TEST.1 REGION: complement(33..53)"
    assert pg.terminator("A") == "TEST.1 REGION: 21..39"
    assert pg.terminator("B") == "TEST.1 REGION: complement(21..39)"
    assert pg.promoter("B") == "TEST.1 REGION: complement(51..69)"
    assert pg.promoter("C") == "TEST.1 REGION: 51..69"
    for function, name in [(pg.promoter, "A"), (pg.terminator, "C")]:
        with pytest.raises(ValueError, match="no neighboring"):
            function(name)


@pytest.mark.parametrize("interval", [[[20, 30]], [[21, 30]], [[5, 25]], [[10, 20]]])
def test_no_gap(synthetic, interval):
    synthetic["genes"]["B"]["gene"] = interval
    with pytest.raises(ValueError, match="overlapping|touching"):
        pg.terminator("A")


def test_nested_neighbor_blocks_its_boundary(synthetic):
    synthetic["genes"]["A"]["gene"] = [[1, 60]]
    with pytest.raises(ValueError, match="overlapping"):
        pg.promoter("B")
    with pytest.raises(ValueError, match="overlapping"):
        pg.terminator("B")


def test_missing_ambiguous_cds_and_alias(synthetic):
    synthetic["genes"]["A"]["cds"] = []
    with pytest.raises(ValueError, match="found 0"):
        pg.cds("A")
    synthetic["genes"]["A"]["cds"] = [[[10, 20]], [[12, 20]]]
    with pytest.raises(ValueError, match="found 2"):
        pg.cds("A")
    with pytest.raises(ValueError, match="A, C"):
        pg.gene("amb")


@pytest.mark.parametrize("function", [pg.gene, pg.cds, pg.locus, pg.promoter, pg.terminator])
def test_invalid_lookup(function):
    for name, error in [(None, TypeError), ("  ", ValueError), ("not_a_gene", KeyError)]:
        with pytest.raises(error):
            function(name)
    with pytest.raises(ValueError, match="Unsupported genome"):
        function("ACT1", "other")
    with pytest.raises(TypeError):
        function("ACT1", None)


@pytest.mark.parametrize("value,error", [(-1, ValueError), (True, TypeError), (1.5, TypeError), ("2", TypeError)])
@pytest.mark.parametrize("parameter", ["upstream", "downstream"])
def test_invalid_flanks(value, error, parameter):
    with pytest.raises(error):
        pg.locus("ACT1", **{parameter: value})


@pytest.mark.parametrize("expression,location", [
    ("1..10", "1:10:1"), ("complement(1..10)", "1:10:2"),
    ("join(1..10,20..30)", "1:10:1,20:30:1"),
    ("complement(join(1..10,20..30))", "20:30:2,1:10:2")])
def test_links(expression, location):
    assert pg.genbanklink("TEST.1 REGION: " + expression) == "https://www.ncbi.nlm.nih.gov/nuccore/TEST.1?location=" + location


def test_act1_link():
    assert pg.genbanklink(pg.cds("ACT1")) == "https://www.ncbi.nlm.nih.gov/nuccore/BK006940.2?location=54687:54696:2,53260:54377:2"


@pytest.mark.parametrize("expression", ["0..10", "10..1", "-1..3", "join(1..2)", "join(5..8,1..2)",
    "join(1..4,4..8)", "complement(join(1..2,4..5)", "1..2 garbage", "<1..2", "order(1..2,3..4)",
    "complement(complement(1..2))", "1..2,3..4", ""])
def test_invalid_links(expression):
    with pytest.raises(ValueError):
        pg.genbanklink("TEST.1 REGION: " + expression)


def test_link_types():
    with pytest.raises(TypeError):
        pg.genbanklink(None)
    with pytest.raises(ValueError):
        pg.genbanklink("https://bad REGION: 1..2")
