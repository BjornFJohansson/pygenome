# This is a plan for a new version of Pygenome 

The new version does not distribute the sequence files of S. cerevisiae, but rather the accession number of the chromosome where
each gene is located and the coordinates.

The python code below shows how it should work.

```python
from pygenome import gene, cds, locus, promoter, terminator, genbanklink

gene("ACT1", genome="S288C") == "BK006940.2 REGION: complement(53260..54696)"  #  has intron
cds("ACT1", genome="S288C") == "BK006940.2 REGION: complement(join(53260..54377,54687..54696))" #  no intron
locus("ACT1", genome="S288C") == "BK006940.2 REGION: complement(52260..55696)" # 1000 up and down string, with intron

genbanklink("BK006940.2 REGION: complement(join(53260..54377,54687..54696))") # returns a string https://www.ncbi.nlm.nih.gov/nuccore/BK006940.2?location=54687:54696:2,53260:54377:2

assert promoter("ACT1") == terminator("YPT1") == "BK006940.2 REGION: complement(54697..55365)"  # strings
```


I would like to have a "scripts" folder containing scripts useful for creating lists and dicts conatining the information from 
the sixteen chromosome files.

For example, the dict sysname:

The dict sysname has a systematic name for each standard name and synonym.

loop over files chr01.gbf .. chr16.gbf in the data folder.

for each file loop over features of type gene:

     gene            complement(53260..54696)
                     /gene="ACT1"
                     /locus_tag="YFL039C"
                     /gene_synonym="ABY1"
                     /gene_synonym="END7"

     mRNA            complement(join(53260..54377,54687..54696))
                     /gene="ACT1"
                     /locus_tag="YFL039C"
                     /gene_synonym="ABY1"
                     /gene_synonym="END7"
                     /product="actin"

Add gene and gene_synonym as keys and locus_tag as values:

for example:  sysname = {"ACT1": "YFL039C", "ABY1": "YFL039C", "END7": "YFL039C" ... }

In the same loop, make the dicts gene and cds:

for example:

genedict = {"YFL039C": (53260..54696), ... }
cdsdict = {"YFL039C": ((53260, 54377), (54687, 54696)), ... }



function:

def gene("ACT1", genome="S288C") -> str: 

    # "BK006940.2 REGION: complement(53260..54696)"  has intron

return acc_region

def cds("ACT1", genome="S288C") == 

    #"BK006940.2 REGION: complement(join(53260..54377,54687..54696))" has  no intron

return acc_region

def locus("ACT1", genome="S288C") == 

    # "BK006940.2 REGION: complement(52260..55696)" # 1000 up and 1000 bp down with intron

return acc_region




Hw should these functions work? 


for example gene("ACT1"):

1. look up sysname

sysname = ["ACT1": "YFL039C", "ABY1": "YFL039C", "END7": "YFL039C" ... ]

2. Look up accession number by the second letter of the systematic name.

accnum = {"F": "BK006940.2", ...}

3. Look up positions

genedict["YFL039C"] == (53260, 54696) 

4. Deduce complement or not from the last letter in sysname, "W" or "C"

return "BK006940.2 REGION: complement(53260..54696)"

---

## Suggested implementation prompt

The following is a proposed specification for implementing the ideas above.
These defaults resolve ambiguities so implementation can proceed independently;
they are suggestions, not claims that the existing examples have been verified.

### Objective and scope

Implement a new coordinate-only Pygenome API in `src/pygenome`, using the existing
`data/chr01.gbf` through `data/chr16.gbf` as the authoritative input. Support
the sixteen nuclear chromosomes of S. cerevisiae S288C initially. Do not bundle
chromosome sequences or download data during import or normal API calls.

Inspect the repository and applicable instructions before editing. This branch
is dedicated to the new implementation. Replace or remove legacy code as needed;
there is no requirement to preserve the old sequence-object API or maintain
backward compatibility. Do not create backup copies or retain old versions of
code. Remove files, tests, documentation, and dependencies that are no longer
needed for the new implementation, after checking their usage. Keep the source
GenBank files and anything required to generate, test, document, or distribute
the new API. This cleanup is explicitly authorized as part of implementation
and does not require separate confirmation.
Implement, document, and verify the result rather than stopping at a proposal.

### Public API

Export these functions directly from `pygenome`:

```python
gene(name: str, genome: str = "S288C") -> str
cds(name: str, genome: str = "S288C") -> str
locus(name: str, genome: str = "S288C", *, upstream: int = 1000,
      downstream: int = 1000) -> str
genbanklink(region: str) -> str
```

Accept systematic names, standard names, and unambiguous synonyms. Strip leading
and trailing whitespace and use case-insensitive lookup. Include systematic
names themselves in the lookup. Raise `KeyError` for an unknown name and
`ValueError` for an ambiguous alias, unsupported genome, or invalid value;
use `TypeError` for incorrect argument types. Reject negative flank lengths,
booleans used as lengths, and empty names. Error messages should identify the
problem and, for ambiguous aliases, list the candidate systematic names.

`gene` returns the annotated gene location, including introns. `cds` returns
the annotated CDS location, excluding introns; use CDS features, not mRNA
features, because transcripts may also contain untranslated regions. Raise
`ValueError` when a gene has no CDS or has multiple distinct CDS locations that
cannot be selected unambiguously. Do not silently pick the first feature.

`locus` extends the outer gene boundaries in the biological upstream and
downstream directions. On the reverse strand, upstream extends toward larger
coordinates. Clip flanks to chromosome boundaries. Include the entire gene,
including introns, and preserve its strand.

### Coordinates and generated data

Use 1-based inclusive coordinates for stored intervals and returned strings.
Convert explicitly from the parser's coordinate convention. Format intervals
as `start..end`, multiple intervals as `join(...)`, and reverse-strand locations
as `complement(...)`. Store joined intervals in ascending genomic order and put
one complement around the whole expression for a reverse-strand feature.

Read accession and version from each GenBank record. Preserve the version in
every result. Store chromosome lengths, strand, and accession alongside each
gene's coordinates. Do not infer chromosome or strand solely from systematic
names: validate the name-based convention against annotations where applicable,
but use the annotations as the source of truth, including for unusual names.

Create `scripts/build_genome_data.py` with explicit input and output path options.
It should generate a deterministic, readable metadata file packaged under
`src/pygenome`, containing name mappings, gene/CDS intervals, chromosome metadata,
and provenance (input filenames, accession versions, and SHA-256 checksums).
Use a single canonical representation rather than duplicating data unnecessarily.
Load packaged resources independently of the working directory.

Parse all sixteen input files. Match gene and CDS features by `locus_tag` and
collect all standard names and repeated synonym qualifiers. Retain alias
collisions as ambiguities; never overwrite them silently. Report missing tags,
unmatched CDS features, duplicate/conflicting records, and unsupported locations.
Do not coerce fuzzy coordinates or unsupported location operators to exact
coordinates. Fail clearly on unrepresentable locations rather than silently
omitting genes. Missing CDS annotations are valid for noncoding genes.

Keep data generation offline and reproducible. A repeated run on the same inputs
must produce byte-identical output. Use Biopython for generation if appropriate,
declaring it explicitly as a development dependency. Prefer the standard library
for the runtime API; remove obsolete dependencies after checking usage.

### GenBank links

Parse and validate the complete region string. Support single intervals,
`join(...)`, and `complement(...)` around either form. Reject malformed syntax,
nonpositive coordinates, and reversed interval bounds. For the reverse strand,
emit joined URL segments in biological order, using strand code `2`; use `1`
for the forward strand. Verify the URL convention before documenting it as
supported, without making runtime calls to NCBI.

The intended corrected example is:

```python
genbanklink("BK006940.2 REGION: complement(join(53260..54377,54687..54696))")
# https://www.ncbi.nlm.nih.gov/nuccore/BK006940.2?location=54687:54696:2,53260:54377:2
```

### Promoters and terminators

Implement and export `promoter(name: str, genome: str = "S288C") -> str` and
`terminator(name: str, genome: str = "S288C") -> str` in this version. These names
refer only to intergenic intervals, not experimentally defined regulatory regions.
Like `gene`, they return accession/region strings rather than sequence objects.

`promoter` returns the interval between the requested gene's upstream boundary
and the nearest neighboring gene in that direction. `terminator` returns the
corresponding interval on its downstream side. Determine upstream and downstream
from the requested gene's strand, regardless of the neighboring gene's strand.
Use the outer boundaries of annotated gene features, not CDS features.
Only genes with annotated CDS features count as neighbors. Ignore noncoding
RNA annotations when selecting neighbors or checking for overlap, including
annotations spanning the gap. Keep noncoding genes available to ordinary lookup.
In particular, `promoter("GAL1")` must return
`"BK006936.2 REGION: 278353..279020"`, and `promoter("GAL10")` must return
`"BK006936.2 REGION: complement(278353..279020)"`, despite the overlapping
noncoding RNA annotation YNCB0008W.

Exclude both genes' bases from the interval. For adjacent genes occupying
`a..b` and `c..d`, where `b < c`, the intergenic interval is `b+1..c-1`.
Return the interval in the requested gene's orientation, using `complement(...)`
for a reverse-strand gene, and preserve the accession version. Apply the same
name resolution and genome validation as the other functions.

If the genes overlap or touch, raise `ValueError` because there is no nonempty
intergenic interval. If no neighboring gene exists before the chromosome end,
raise `ValueError` because the interval is not bounded by two genes. Do not skip
an overlapping neighbor to select a more distant gene.

Verify the opening `promoter("ACT1") == terminator("YPT1")` example against the
source annotations. If these genes bound the same gap and have the same strand,
the strings should be identical. Correct the example's coordinates to exclude
gene bases and include the accession version. Add tests for both strands,
neighbors on opposite strands, overlaps, touching genes, and chromosome ends.

### Verification and delivery

Add focused tests for forward/reverse strands, intron-containing and intronless
genes, noncoding genes, synonyms and collisions, unequal flanks, chromosome
boundaries, missing/ambiguous CDS annotations, and malformed region strings.
Use small synthetic GenBank fixtures for edge cases and real input files for
integration checks. Verify the ACT1 examples against the annotations; if an
example conflicts with the input data, report the discrepancy and use the data.

Check every generated interval against chromosome bounds, account for all input
gene features, and verify deterministic regeneration. Build a wheel and source
distribution, check that metadata is included and chromosome sequences are
excluded, and smoke-test the installed API from outside the checkout without
network access. Update or remove obsolete legacy tests so the maintained test
suite validates the new API and passes.

Write a README with installation, API examples, coordinate conventions, errors,
data provenance, and the regeneration command. Finish with a concise report of
changed files, checks run, and any unresolved annotation issues. Do not commit
or publish unless requested.

Proceed using these defaults without asking about routine implementation choices.
Ask only if essential source data is unavailable or a biological ambiguity would
require changing this specification. Otherwise document the decision and continue.
