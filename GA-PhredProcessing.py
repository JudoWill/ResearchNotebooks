# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import os, os.path
import shutil
import glob
import sys
from subprocess import check_call, check_output

os.chdir('/home/will/Dropbox/PhredDirectory/')
staden_path = '/home/will/staden-2.0.0b9.x86_64/bin/'
sys.path.append('/home/will/PySeqUtils/')

# <codecell>

from GeneralSeqTools import call_muscle, fasta_reader, fasta_writer

# <codecell>

#from Bio import SeqIO
from Bio.SeqIO.AbiIO import AbiIterator
files = glob.glob('../Wigdahl Trace files/2:11:11/*.ab1')
seqs = []
for f in files:
    rec = AbiIterator(open(f, mode = 'rb'), trim = True).next()
    seqs.append( (rec.id, rec.seq.tostring()) )

# <codecell>

!/home/will/staden-2.0.0b9.x86_64/bin/convert_trace --help

# <codecell>

res = call_muscle(seqs)
with open('align_data.fasta', 'w') as handle:
    fasta_writer(handle, res)

# <codecell>

from HIVTransTool import process_seqs

results = list(process_seqs(seqs[:50], extract_regions = True, known_names = 50))

# <codecell>

for row in results:
    if row['RegionName'] == 'LTR5':
        print row['Name'], row['QueryNuc']

# <codecell>

results[:5]

# <codecell>


