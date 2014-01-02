# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import sys
sys.path.append('/home/will/PySeqUtils/')
from GeneralSeqTools import fasta_reader, fasta_writer
import os
os.chdir('/home/will/PySeqUtils/TransToolStuff/')

# <codecell>

from itertools import islice
start = 806
stop = -1
path = 'HIV1_ALL_2012_env_PRO.fasta'
outpath = 'HIV1_ALL_2012_gp41_PRO.fasta'
with open(path) as handle:
    for name, seq in islice(fasta_reader(handle), 20):
        tseq = seq[start:stop] 
        print tseq[:5], tseq[-5:]

# <codecell>

seqs = []
with open(path) as handle:
    for name, seq in fasta_reader(handle):
        seqs.append((name, seq[start:stop]))
with open(outpath, 'w') as handle:
    fasta_writer(handle, seqs)

# <codecell>

from Bio import Entrez
from Bio import SeqIO
ids = '544451412,544451410,544451408,544451406,544451404,544451402,544451400,544451398,544451396'

fetch_handle = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text",
                             id=ids)
records = list(SeqIO.parse(fetch_handle, "gb"))

# <codecell>

rec = records[0]

# <codecell>

rec.annotations['gi']

# <codecell>

batch_size = 1000
num_res = 1300
inds = range(0, num_res, batch_size)+[num_res]
start_inds = inds
stop_inds = inds[1:]
zip(start_inds, stop_inds)

# <codecell>

from collections import defaultdict

counts = defaultdict(int)
seqs = []
with open('/home/will/WLAHDB_data/RegionDBs/LTR/HIV1_ALL_2012_ltr_DNA.fasta') as handle:
    for name, seq in fasta_reader(handle):
        seqs.append((name, seq.replace('-', '')))

with open('/home/will/WLAHDB_data/RegionDBs/LTR/LTR.fasta', 'w') as handle:
    fasta_writer(handle, seqs)

# <codecell>

max(len(s) for _, s in seqs)

# <codecell>

from Bio import SeqIO
from Bio.Alphabet import generic_dna
with open('/home/will/WLAHDB_data/SubtypeDB/HIV1_genome_DNA.fasta') as handle:
    seqs = list(SeqIO.parse(handle, 'fasta', generic_dna))
    

# <codecell>

import glob

files = glob.glob('/home/will/WLAHDB_data/GenbankDL/*.gb')

counts = []
for fnum, f in enumerate(files):
    if fnum % 10000 == 0:
        print fnum, sum(counts)
    with open(f) as handle:
        for num, seq in enumerate(SeqIO.parse(handle, 'gb'), 1):
            pass
        counts.append(num)


# <codecell>

from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastxCommandline, NcbiblastnCommandline
from tempfile import NamedTemporaryFile
from StringIO import StringIO

window_size = 500
inds = range(0, len(seq), window_size) + [len(seq)]
seq = seqs[0]
blocks = []
for start, stop in zip(inds, inds[1:]):
    blocks.append(seq[start:stop])

with NamedTemporaryFile(suffix='.fasta', delete=False) as handle:
    with NamedTemporaryFile() as ohandle:
    
        SeqIO.write(blocks, handle, 'fasta')
        handle.flush()
        os.fsync(handle.fileno())
    
        cline = NcbiblastnCommandline(db='/home/will/WLAHDB_data/SubtypeDB/HIV1_genome_DNA.fasta',
                                      query=handle.name,
                                      out=ohandle.name,
                                      outfmt=5)
        _, _ = cline()
        records = list(NCBIXML.parse(ohandle))


# <codecell>

rec = records[0]
align = rec.alignments[0]

# <codecell>

rec.query

# <codecell>

align.hit_def

# <codecell>


