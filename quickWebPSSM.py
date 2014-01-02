# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import sys
import os, os.path
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
from itertools import chain
from pandas import DataFrame, Series, merge, read_csv, MultiIndex, Index, concat

sys.path.append('/home/will/PySeqUtils/')
from GeneralSeqTools import fasta_reader, WebPSSM_V3_fasta, yield_chunks

# <codecell>

seq_path = '/home/will/HIVTropism/Harrigan_2012_V3.fasta'
with open(seq_path) as handle:
    v3_seqs = list(fasta_reader(handle))

have_data = []
need_data = set([name for name, _ in v3_seqs])
count = 0
print 'Getting WebPSSM scores'
while need_data and count < 5:
    count += 1
    print len(need_data)
    gen = ((name, seq) for name, seq in v3_seqs if name in need_data)
    chunks = yield_chunks(gen, 50)
    with ThreadPoolExecutor(max_workers = 10) as e:
        res = e.map(WebPSSM_V3_fasta, chunks)
        have_data += list(chain.from_iterable(res))
    need_data -= set(row[0] for row in have_data)

# <codecell>

fields = ['name','score','pred','x4.pct','r5.pct','geno','pos.chg','net.chg','percentile']
pssm_data = DataFrame(have_data, columns=fields)

# <codecell>

outpath = '/home/will/HIVTropism/Harrigan_2012_V3.csv'
pssm_data.to_csv(outpath, index=False)

