"""S288C gene coordinates (1-based, inclusive)."""
from functools import lru_cache
from importlib.resources import files
import json
import re

__all__ = ["gene", "cds", "locus", "promoter", "terminator", "genbanklink"]


@lru_cache(maxsize=1)
def _data():
    return json.loads(files(__package__).joinpath("S288C.json").read_text(encoding="utf-8"))


def _resolve(name, genome):
    if not isinstance(genome, str):
        raise TypeError("genome must be a string")
    if genome != "S288C":
        raise ValueError(f"Unsupported genome: {genome!r}; expected 'S288C'")
    if not isinstance(name, str):
        raise TypeError("name must be a string")
    key = name.strip().upper()
    if not key:
        raise ValueError("name must not be empty")
    data = _data()
    candidates = data["aliases"].get(key)
    if not candidates:
        raise KeyError(f"Unknown gene: {name!r}")
    if len(candidates) != 1:
        raise ValueError(f"Ambiguous alias {name!r}: {', '.join(candidates)}")
    tag = candidates[0]
    return tag, data["genes"][tag]


def _format(entry, intervals):
    expression = ",".join(f"{a}..{b}" for a, b in intervals)
    if len(intervals) > 1:
        expression = f"join({expression})"
    if entry["strand"] == -1:
        expression = f"complement({expression})"
    return f"{entry['accession']} REGION: {expression}"


def gene(name: str, genome: str = "S288C") -> str:
    """Return the annotated gene location, including introns."""
    _, entry = _resolve(name, genome)
    return _format(entry, entry["gene"])


def cds(name: str, genome: str = "S288C") -> str:
    """Return the unique annotated CDS; reject missing or ambiguous CDSs."""
    tag, entry = _resolve(name, genome)
    if len(entry["cds"]) != 1:
        raise ValueError(f"{tag}: expected one CDS location, found {len(entry['cds'])}")
    return _format(entry, entry["cds"][0])


def locus(name: str, genome: str = "S288C", *, upstream: int = 1000,
          downstream: int = 1000) -> str:
    """Extend the gene in biological directions, clipping to chromosome bounds."""
    _, entry = _resolve(name, genome)
    for label, value in (("upstream", upstream), ("downstream", downstream)):
        if type(value) is not int:
            raise TypeError(f"{label} must be an integer (not a boolean)")
        if value < 0:
            raise ValueError(f"{label} must be nonnegative")
    left, right = (upstream, downstream) if entry["strand"] == 1 else (downstream, upstream)
    length = _data()["chromosomes"][entry["accession"]]["length"]
    return _format(entry, [[max(1, entry["gene"][0][0] - left),
                             min(length, entry["gene"][-1][1] + right)]])


def _intergenic(name, genome, upstream):
    tag, entry = _resolve(name, genome)
    start, end = entry["gene"][0][0], entry["gene"][-1][1]
    left = (entry["strand"] == 1) == upstream
    neighbors = [v for k, v in _data()["genes"].items()
                 if k != tag and v["accession"] == entry["accession"] and v["cds"]]
    # An annotation crossing the queried boundary precludes an intergenic gap.
    boundary = start if left else end
    if any(v["gene"][0][0] <= boundary <= v["gene"][-1][1] for v in neighbors):
        raise ValueError(f"{tag}: overlapping gene at requested boundary")
    if left:
        edges = [v["gene"][-1][1] for v in neighbors if v["gene"][-1][1] < start]
        interval = [max(edges) + 1, start - 1] if edges else None
    else:
        edges = [v["gene"][0][0] for v in neighbors if v["gene"][0][0] > end]
        interval = [end + 1, min(edges) - 1] if edges else None
    if interval is None:
        raise ValueError(f"{tag}: no neighboring gene before chromosome end")
    if interval[0] > interval[1]:
        raise ValueError(f"{tag}: touching genes leave no intergenic interval")
    return _format(entry, [interval])


def promoter(name: str, genome: str = "S288C") -> str:
    """Return the gap to the upstream coding gene, excluding gene bases."""
    return _intergenic(name, genome, True)


def terminator(name: str, genome: str = "S288C") -> str:
    """Return the gap to the downstream coding gene, excluding gene bases."""
    return _intergenic(name, genome, False)


def genbanklink(region: str) -> str:
    """Convert an exact accession/region string to an NCBI nuccore link."""
    if not isinstance(region, str):
        raise TypeError("region must be a string")
    match = re.fullmatch(r"([A-Za-z0-9_]+(?:\.\d+)?) REGION: (.+)", region)
    if not match:
        raise ValueError(f"Malformed region: {region!r}")
    accession, expression = match.groups()
    reverse = expression.startswith("complement(") and expression.endswith(")")
    if reverse:
        expression = expression[11:-1]
    joined = expression.startswith("join(") and expression.endswith(")")
    if joined:
        expression = expression[5:-1]
    pattern = r"[1-9]\d*\.\.[1-9]\d*"
    if not re.fullmatch(pattern + (rf"(?:,{pattern})+" if joined else ""), expression):
        raise ValueError(f"Unsupported or malformed location: {region!r}")
    intervals = [tuple(map(int, piece.split(".."))) for piece in expression.split(",")]
    if any(a > b for a, b in intervals) or any(a[1] >= b[0] for a, b in zip(intervals, intervals[1:])):
        raise ValueError("Intervals must be ascending, disjoint, and have start <= end")
    if reverse:
        intervals.reverse()
    location = ",".join(f"{a}:{b}:{2 if reverse else 1}" for a, b in intervals)
    return f"https://www.ncbi.nlm.nih.gov/nuccore/{accession}?location={location}"
