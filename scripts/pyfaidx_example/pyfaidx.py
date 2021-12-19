from pyfaidx import Fasta

chromosome = Fasta('chr01.fsa')

chromosome.keys()

x = chromosome['BK006935.2'][5:10]

x.complement


y = -x
