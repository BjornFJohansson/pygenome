#!/usr/bin/env python3

"""
gene("ACT1")

1. look up sysname

sysname = ["ACT1": "YFL039C", "ABY1": "YFL039C", "END7": "YFL039C" ... ]

2. Look up accession number by the second letter of the systematic name.

accnum = ["F": "BK006940.2"]

3. Look up positions

gene_position["YFL039C": (53260, 54696),  ... ]

deduce complement or not from the last letter in sysname.

return "BK006940.2 REGION: complement(53260..54696)"

cds("ACT1") = 

Deletion_primers_PCR_sizes.txt: line 1885, second column = "YFL039C"

"""

# YFL039C

from pygenome import gene, cds, locus, genbanklink

gene("ACT1", genome="S288C") == "BK006940.2 REGION: complement(53260..54696)"  # string, has intron
cds("ACT1", genome="S288C") == "BK006940.2 REGION: complement(join(53260..54377,54687..54696))" # string, no intron
locus("ACT1", genome="S288C") == "BK006940.2 REGION: complement(52260..55696)" # 1000 up and down (string), with intron

genbanklink("BK006940.2 REGION: complement(join(53260..54377,54687..54696)") # returns a string https://www.ncbi.nlm.nih.gov/nuccore/BK006940.2?location=54687:54696:2,53260:54377:2

ddbjlink("....")   # returns a string
enalink("....")    # returns a string
sgdlink("....")    # returns a string

from pydna.genbank import genbank
act1orf = genbank("BK006940.2 REGION: complement(53260..54696)", email="bjornjobb@gmail.com")  # Dseqrecord
act1locus = genbank("BK006940.2 REGION: complement(52260..55696)", email="bjornjobb@gmail.com")  # Dseqrecord

assert promoter("ACT1") == terminator("YPT1") == "BK006940 REGION: complement(54696..55366)"  # strings

from pydna.genbank import genbank

genbank("BK006940.2 REGION: complement(53260..54696)", email="bjornjobb@gmail.com")

sf = accession_region_to_sf("BK006940.2 REGION: complement(53260..54696)") # biopython sequence feature

from pydna-utils import local_genbank

local_genbank.fetch.url("http://downloads.yeastgenome.org/sequence/S288C_reference/NCBI_genome_source/chr06.gbf")  # fetch a genbank file from genbank
local_genbank.fetch.genbank("BK006940.2")  # fetch a genbank file from genbank
local_genbank.list()  # list the genbank files and their sizes and ACCESSSIONs
dsr = local_genbank("BK006940.2")

sf.extract(dsr)


# https://www.yeastgenome.org/locus/ACT1/sequence
```


