# Pygenome

Offline accession and coordinate lookup for the sixteen nuclear chromosomes of
*S. cerevisiae* S288C. Returns GenBank region strings compatible with pydna and 
pydna-utils for quicker access by caching downloaded sequences.


Requires Python 3.12.7 or newer. Install from this checkout with `pip install .`.
For development, run `uv sync --frozen`, then `uv run pytest`.

```python
from pygenome import gene, cds, locus, promoter, terminator, genbanklink

assert gene("ACT1") == "BK006940.2 REGION: complement(53260..54696)"
assert cds("ACT1") == "BK006940.2 REGION: complement(join(53260..54377,54687..54696))"
assert locus("ACT1") == "BK006940.2 REGION: complement(52260..55696)"
assert promoter("ACT1") == terminator("YPT1") == "BK006940.2 REGION: complement(54697..55365)"
assert genbanklink(cds("ACT1")) == "https://www.ncbi.nlm.nih.gov/nuccore/BK006940.2?location=54687:54696:2,53260:54377:2"
```

All lookup functions accept `name` and optional `genome="S288C"`. Names may be
systematic names, standard names, or synonyms; lookup ignores case
and surrounding whitespace. Only the exact genome identifier `S288C` is supported.

| Function | Result |
| --- | --- |
| `gene(name, genome="S288C")` | Annotated gene location, including introns |
| `cds(name, genome="S288C")` | Annotated CDS location, excluding introns |
| `locus(name, genome="S288C", *, upstream=1000, downstream=1000)` | Outer gene span plus strand-aware flanks, clipped at chromosome ends |
| `promoter(name, genome="S288C")` | Upstream intergenic interval |
| `terminator(name, genome="S288C")` | Downstream intergenic interval |
| `genbanklink(region)` | NCBI nuccore URL for an exact region string |

Coordinates are **1-based and inclusive**. Joined intervals are ordered by
ascending genomic position, with one outer `complement(...)` on the reverse
strand. Accession versions are preserved. Chromosome and strand come from the
annotations, rather than being inferred from names.

Promoter and terminator mean only the gaps between annotated genes. They exclude
both genes' bases and use the requested gene's orientation, regardless of the
neighbor's strand. Only genes with annotated CDS features count as neighbors;
use their full gene boundaries, not CDS boundaries. Noncoding RNA annotations
do not interrupt these intervals, even when they overlap a gene or span the gap.

For example, GAL1 and GAL10 share `278353..279020` on `BK006936.2`:
`promoter("GAL1")` returns the forward interval and `promoter("GAL10")` its
complement. Noncoding genes remain available to all lookup functions, but their
promoter/terminator neighbors are also selected from coding genes.
A coding neighbor crossing the requested gene boundary blocks that side;
a gene wholly contained within the requested gene does not alter its external
gaps. Missing neighbors at chromosome ends, touching genes, and overlapping
boundaries produce errors rather than empty or invented intervals.

Unknown names raise `KeyError`. Ambiguous aliases (with candidate names in the
message), unsupported genomes, empty names, negative flanks, missing/multiple
CDS locations, unavailable intergenic intervals, and malformed region strings
raise `ValueError`. Incorrect types, including boolean flank lengths, raise
`TypeError`. The current dataset contains 160 ambiguous aliases.

`genbanklink` accepts exact intervals, joins, and complements of either. Joined
intervals must be ascending and disjoint. It emits reverse-strand segments in
biological order with strand code `2` (`1` for the forward strand), following
the URL convention specified in `plan.md`. 

Live verification of joined-region rendering on NCBI was blocked by NCBI's browser 
challenge during implementation; the URL construction is tested locally.

## Data and regeneration

The authoritative inputs are `data/chr01.gbf` through `data/chr16.gbf`. Their
source URLs are retained in `data/chromosome_urls.txt`. The generated
`src/pygenome/S288C.json` records each filename, versioned accession, SHA-256
checksum, chromosome length, gene strand and intervals, CDS intervals, aliases,
and an annotation audit. It contains no nucleotide or protein sequences.

The current snapshot has 6,424 genes and 6,001 CDS features; 423 genes have no
CDS. Every input gene is accounted for. There are no systematic-name chromosome
or strand mismatches and no multiple distinct CDS locations in this snapshot.

From a checkout containing the source GenBank files:

```sh
uv sync --frozen
uv run python scripts/build_genome_data.py --input-dir data --output src/pygenome/S288C.json
uv run pytest
uv build
```

Generation is offline and deterministic. Biopython is a development dependency,
not a runtime dependency. The generator checks all gene and CDS annotations,
rejects fuzzy/unsupported locations, missing tags, duplicate genes and unmatched
or inconsistent CDS features, and retains alias collisions explicitly. Repeated
identical CDS locations are deduplicated; distinct locations are retained so
`cds` can report ambiguity. Generated chromosome lengths are normalized in the
chromosome table, referenced by each gene's accession.

The source distribution includes the generator and synthetic tests but excludes
the large GenBank files. Its source-data integration test skips when those inputs
are absent. To regenerate metadata, use the repository checkout with its matching
input files. The legacy sequence-object API and deletion-cassette workflows have
been removed for this implementation.
