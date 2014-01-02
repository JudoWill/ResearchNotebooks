# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import os, os.path
import shutil
import glob

# <codecell>

os.path.join('home', 'will', 'Dropbox', 'IpythonShare')

# <codecell>

#list all files in a path
flist = os.listdir(os.path.join('/home', 'will', 'Dropbox'))
glob.glob(os.path.join('/home', 'will', 'Dropbox', '*','*.txt'))

# <codecell>

import csv

fname = '/home/will/HIVTropism/more_phylip_BenjRes.tsv'
with open(fname) as handle:
    reader = csv.DictReader(handle, delimiter='\t')
    count = 0
    for row in reader:
        print row
        count += 1
        if count > 10:
            break

# <codecell>


