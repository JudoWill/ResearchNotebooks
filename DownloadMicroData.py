# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import os, os.path
import concurrent.futures
import csv
import urllib.request
import shutil
import gzip

os.chdir('/home/will/AutoMicroAnal/')

# <codecell>

with open('MicroarraySamples.tsv') as handle:
    microdata = list(csv.DictReader(handle, delimiter = '\t'))

# <codecell>

def get_fileurl(supurl):
    #print(supurl)
    resp = urllib.request.urlopen(supurl)
    for line in resp:
        fname = str(line.split()[-1])
        if fname.lower().endswith(".cel.gz'"):
            #print('returning')
            return supurl + fname[2:-1]
    return None

def process_row(row):
    supurl = row['URL'] + 'suppl/'
    tmpfile = '/tmp/' + row['Sample Accession'] + '.CEL.gz'
    finalfile = '/home/will/AutoMicroAnal/microadata/' + row['Sample Accession'] + '.CEL'
    if os.path.exists(finalfile):
        return None
    fileurl = get_fileurl(supurl)
    #print(fileurl)
    if fileurl is None:
        return fileurl
    try:
        resp = urllib.request.urlopen(fileurl)
        with open(tmpfile, 'wb') as handle:
            handle.write(resp.read())
    except urllib.request.URLError:
        return fileurl
        
    with gzip.open(tmpfile) as zhandle:
        with open(finalfile, 'wb') as handle:
            handle.write(zhandle.read())
    os.remove(tmpfile)
    return None

# <codecell>

gp133As = [row for row in microdata if row['Platform'] == 'GPL96']

# <codecell>

for num, row in enumerate(gp133As):
    try:
        res = process_row(row)
    except:
        print('skipping')
        continue
    if (num == 0) | (num == 5) | (num == 20) | (num % 500 == 0):
        print(num)
    if res:
        print(res)

# <codecell>


# <codecell>


