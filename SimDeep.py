# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import os, os.path
import random
os.chdir('/home/will/DeepSequencing/testfix/')

# <codecell>

with open('RefSequence_single.fasta') as handle:
    _ = handle.next()
    seq = handle.next().strip()[:600]

# <codecell>

haplotypes = []
for _ in range(100):
    tmp = range(len(seq))
    random.shuffle(tmp)
    nseq = list(seq)
    for n in tmp[:50]:
        nseq[n] = 'A'
    haplotypes.append(''.join(nseq))

# <codecell>

tfracs = [0.4, 0.2, 0.1, 0.05]
nleft = len(haplotypes) - len(tfracs)
fracs = tfracs + [(1-sum(tfracs))/(nleft)]*nleft
cfracs = [sum(fracs[:n]) for n in range(1,len(fracs))]
sum(fracs)
print cfracs

# <codecell>

seq_len = 80
nreads = 500000
pos_sites = range(600-80)
with open('sim_reads.fa', 'w') as handle:
    for _ in range(nreads):
        val = random.uniform(0,1)
        for hap, cval in enumerate(cfracs):
            if val <= cval:
                break
        pos = random.randint(0, 600-80)
        read = haplotypes[hap][pos:pos+seq_len]
        rname = 'hap-%i-start-%i' % (hap, pos)
        handle.write('>%s\n%s\n' % (rname, read))
    
with open('known_haps.fasta', 'w') as handle:
    for hap, (seq, freq) in enumerate(zip(haplotypes, fracs)):
        handle.write('>hap-%i-%f\n%s\n' % (hap, freq, seq))

# <codecell>


