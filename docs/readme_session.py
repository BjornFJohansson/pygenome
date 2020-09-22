# In[]
from pygenome import sg
# In[]
mygene = sg.stdgene["XKS1"]
# In[]
mygene
# In[]
mygene.short_description
# In[]
sg.sysgene["YGR194C"]
# In[]
mygene.cds
# In[]
mygene.locus()
# In[]
mygene.promoter
# In[]
mygene.terminator
# In[]
mygene.downstream_gene
# In[]
mygene.upstream_gene
# In[]
mygene.deletion_locus