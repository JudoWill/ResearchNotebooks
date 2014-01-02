# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import os, os.path
import sys

sys.path.append('/home/will/PySeqUtils/')
os.chdir('/home/will/PySeqUtils/')
from GeneralSeqTools import fasta_reader
import csv

# <codecell>

from HIVTransTool import map_seqs_to_ref, process_seqs
import csv

ref_path = 'HIVDBFiles/HXB2Sequence.fasta'
cor_res = {}
with open('TestData/test_mapping.csv') as handle:
    for row in csv.DictReader(handle, delimiter='\t'):
        cor_res[row['SeqName']] = (int(row['GenomeStart']), int(row['GenomeStop']))

with open('TestData/testSeqs.fasta') as handle:
    input_seqs = list(fasta_reader(handle))
    
#with open('/home/will/HIVRate/hiv-db.fasta') as handle:
#    seqs = list(fasta_reader(handle))

# <codecell>

from itertools import product

list(product('abcdefg', range(5)))

# <codecell>


fields = ['Name','RegionName', 'QueryNucStart','QueryNucStop','QueryNuc',
'RegionNucStart','RegionNucStop','RegionAAStart', 'RegionAAStop', 'QueryAA']
with open('TestData/LocatorRes.tsv', 'w') as handle:
    writer = csv.DictWriter(handle, fields, delimiter='\t')
    writer.writeheader()
    for row in process_seqs(input_seqs, extract_regions = True):
        writer.writerow(row)
        #print row
            
#raise KeyError
    
        

# <codecell>

import csv
headers = ['Name', 'RegionName', 'QueryNucStart', 'QueryNucStop', 'QueryNuc',
           'RegionNucStart', 'RegionNucStop', 'RegionAAStart','RegionAAStop','QueryAA']
with open('test_output.tsv', 'w') as handle:
    writer = csv.DictWriter(handle, headers, delimiter='\t')
    found_pats = set()
    for num, row in enumerate(process_seqs(seqs)):
        found_pats.add(row['Name'])
        if (num == 1) or (num == 100) or (num == 1000) or (num % 5000 == 0):
            print num, len(found_pats)
        writer.writerow(row)



# <codecell>

len(seqs)

# <codecell>


