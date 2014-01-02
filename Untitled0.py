# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import os, os.path
import sys
import pandas as pd
import numpy as np
sys.path.append('/home/will/PatientPicker/')
sys.path.append('/home/will/PySeqUtils/')
import LoadingTools

# <codecell>

redcap_data = LoadingTools.load_redcap_data()

# <codecell>

def get_months(dt):
    
    return dt.astype('timedelta64[D]')/np.timedelta64(1, 'D')/30

date_data = redcap_data[['Patient ID', 'VisitNum','Date Of Visit']]
l = []
for tup in pat_joins:
    
    for pat in tup:
        mask = date_data['Patient ID'] == pat
        tmp = date_data[mask].copy()
        tmp.sort('VisitNum')
        tmp['Date'] = pd.to_datetime(tmp['Date Of Visit'], coerce = True)
        tmp = tmp.drop(['Date Of Visit'], axis = 1)
        try:
            tmp['Months from R00'] = (tmp['Date'] - tmp['Date'].iloc[0]).map(get_months)
        except IndexError:
            pass
        l.append(tmp.copy())
    
    
extra_res = pd.concat(l, axis = 0, ignore_index = True)
extra_res.to_excel('/home/will/Downloads/joining_patients.xlsx')

# <codecell>

pat_joins = [('A0138', 'A0360'),
             ('A0017', 'A0054'),
             ('A0418', 'A0423'),
             ('A0234', 'A0266'),
             ('A0410', 'A0424'),
             ('A0082', 'A0346'),
             ('A0077', 'A0370'),
             ('A0098', 'A0126'),
             ('A0023', 'A0114'),
             ('A0131', 'A0157'),
             ('A0059', 'A0403'),
             ('A0016', 'A0053'),
             ('A0219', 'A0442')]
#pat_index = redcap_data['Patient ID']
#orig = pat_index.copy()
#for joins in pat_joins:
#    base_pat = joins[0]
#    for join_pat in joins[1:]:
#        print base_pat, join_pat, (pat_index == join_pat).sum()
#        pat_index[pat_index == join_pat] = base_pat
#(orig != pat_index).sum()


# <codecell>

import glob
ld = []
for f in glob.glob('/home/will/SubCData/LANLRes/LANLResults_*.txt'):
    print f
    ld.append(pd.read_csv(f, sep='\t', index_col=0))
print 'concating'    
lanl_data = pd.concat(ld, axis = 0, ignore_index=True)

# <codecell>

from types import StringType
def norm_id(inser):
    
    cols = ['Accession', 'GI number']
    for col in cols:
        if type(inser[col]) == StringType:
            return inser[col].split('.')[0]
    print inser
    raise KeyboardInterrupt
    

lanl_data['GBid'] = lanl_data.apply(norm_id, axis = 1)
agg_lanl_data = lanl_data.groupby('GBid').first()

# <codecell>

import glob
from Bio import SeqIO
from concurrent.futures import ProcessPoolExecutor
from itertools import islice, imap
import os, os.path
import csv

def get_gi_acc(fname):
    gb = fname.split('/')[-1].split('.')[0]
    with open(fname) as handle:
        for line in handle:
            if line.startswith('ACCESSION'):
                acc = line.strip().split()[-1]
                return gb, acc
    raise AssertionError



gi_to_acc_dict = {}

fname = '/home/will/WLAHDB_data/gi_to_acc.csv'
if os.path.exists(fname):
    with open(fname) as handle:
        for row in csv.reader(handle):
            gi_to_acc_dict[row[0]] = row[1]
else:
    gb_files = glob.glob('/home/will/WLAHDB_data/GenbankDL/*.gb')
    with open(fname, 'w') as handle:
        writer = csv.writer(handle)
        for num, (gbm, acc) in enumerate(imap(get_gi_acc, gb_files)):
            if (num == 100) or (num % 50000 == 0):
                print num
            gi_to_acc_dict[gbm] = acc
            writer.writerow((gbm, acc))

# <codecell>

import GeneralSeqTools
files = [('C', sorted(glob.glob('/home/will/WLAHDB_data/SeqDump/C_*'))),
         ('B', sorted(glob.glob('/home/will/WLAHDB_data/SeqDump/B_*')))]
seqs = []
for sub, sfiles in files:
    for f in sfiles:
        with open(f) as handle:
            base_name = f.rsplit(os.sep,1)[1].rsplit('.',1)[0]
            prot = base_name.split('_')[1]
            for name, seq in GeneralSeqTools.fasta_reader(handle):
                seqs.append({
                             'Seq':seq,
                             'ID':gi_to_acc_dict[name],
                             'Prot':prot,
                             'Subtype':sub
                             })
            
seqdf = pd.DataFrame(seqs)

# <codecell>

seqdf['HasSeq'] = 1.0
cols = ['CD4 count', 'Viral load']
seq_counts = pd.pivot_table(seqdf,
                            rows = ['Subtype', 'ID'],
                            cols = 'Prot',
                            values = 'HasSeq',
                            aggfunc = 'sum')
ex_data, seq_count_aln = agg_lanl_data[cols].align(seq_counts.ix['C'], axis = 0, join = 'right')

# <codecell>

fourkb_cols = ['vpr', 'tat', 'vpu', 'env', 'ltr']

ex_data['CD4'] = ex_data['CD4 count'].notnull()
ex_data['VL'] = ex_data['Viral load'].notnull()
ex_data['HasV3'] = seq_count_aln['v3'].notnull()
ex_data['Has44'] = seq_count_aln[fourkb_cols].fillna(0).sum(axis = 1)==5

for col in seq_count_aln.columns:
    ex_data[col] = seq_count_aln[col].notnull()

# <codecell>

ex_data[other_cols].sum(axis=1)

# <codecell>

from itertools import product
ncounts = []
other_cols = ['CD4', 'VL', 'HasV3', 'Has44']
for col, ocol in product(seq_count_aln.columns, other_cols+['ASeqs', '_All']):
    if ocol == 'ASeqs':
        ncounts.append({
                    'Prot': col,
                    'Col': ocol,
                    'Count': ex_data[col].sum()
                    })
    elif ocol == '_All':
        ncounts.append({
                    'Prot': col,
                    'Col': ocol,
                    'Count': (ex_data[col]*(ex_data[other_cols].sum(axis=1)==4)).sum()
                    })
    else:
        ncounts.append({
                        'Prot': col,
                        'Col': ocol,
                        'Count': (ex_data[col]*ex_data[ocol]).sum()
                        })
table = pd.pivot_table(pd.DataFrame(ncounts),
               rows = 'Prot',
               cols = 'Col',
               values = 'Count',
               aggfunc = 'sum')
table.to_excel('/home/will/Dropbox/Wigdahl HIV Lab/SubCAnalysis/new_count_table.xlsx')

# <codecell>

table


# <codecell>


# <codecell>

wanted_pvs = [('A0001','R02'),
              ('A0004','R02'),('A0004','R07'),
              ('A0023','R01'),('A0041','R02'),
              ('A0046','R02'),('A0048','R02'),
              ('A0052','R01'),('A0070','R01'),
              ('A0071','R01'),('A0072','R01'),
              ('A0107','R05'),('A0110','R02'),
              ('A0111','R01'),('A0127','R01'),
              ('A0128','R00'),('A0139','R00'),
              ('A0143','R00'),('A0152','R00'),
              ('A0164','R00'),('A0165','R00'),
              ('A0172','R00'),('A0208','R00'),
              ('A0213','R00'),('A0225','R00'),
              ('A0252','R00'),('A0259','R00'),
              ('A0278','R00'),('A0281','R00'),
              ('A0310','R00'),('A0326','R00'),
              ('A0365','R00'),('A0367','R05'),
              ('A0377','R00'),('A0378','R00'),
              ('A0396','R00'),('A0403','R01'),
              ('A0421','R00')]

col_order = ['Patient ID','VisitNum','Age','NumTotalVisits','Latest CD4 count','Current Alcohol Use',
             'Current Tobacco Use','Days since baseline','Gender','Race','HAART','Hepatitis C status (HCV)',
             'HIVD score','HIVD.I','Hepatitis B status (HBV)','Latest CD8 count','Nadir CD4 count',
             'Nadir CD8 count','Peak viral load','Latest Viral load','Years Seropositive',
             'TOSample.Benzodiazepines','TOSample.Cannabinoid','TOSample.Cocaine','TOSample.Opiates',
             'ALL.Benzodiazepines','ALL.Cannabinoid','ALL.Cocaine','ALL.Opiates',
             'ATSample.Amphetamines','ATSample.Barbiturates','ATSample.Benzodiazepines',
             'ATSample.Cannabinoid','ATSample.Cocaine','ATSample.Opiates','ATSample.Phencyclidine']

wanted_pvs = pd.DataFrame(wanted_pvs, columns = ['Patient ID', 'VisitNum'])

# <codecell>

test_cols = [col for col in redcap_data.columns if col.startswith('Test-')]
at_sample = redcap_data[test_cols].rename(dict([(col, col.replace('Test-', 'ATSample.'))]))
to_sample = redcap_data.groupby(level = ['Patient ID'])[test_cols].transform(pd.expanding_mean).rename(dict([(col, col.replace('Test-', 'TOSample.'))]))
all_sample = redcap_data.groupby(level = ['Patient ID'])[test_cols].transform(lambda x: x.mean()).rename(dict([(col, col.replace('Test-', 'ALL.'))]))

# <codecell>

def mark_haart(inrow):
    if inrow['HAART-Naive']:
        return 'nH'
    elif inrow['HAART-Non-Adherent']:
        return 'dH'
    elif inrow['HAART-Off']:
        return 'dH'
    elif inrow['HAART-On']:
        return 'cH'
    else:
        return np.nan
    
    
haart_cols = [col for col in redcap_data.columns if col.startswith('HAART')]
haart_data = redcap_data.groupby(level = [0,1])[haart_cols].apply(mark_haart)
num_data = redcap_data.groupby(level = [0])[['Age']].transform(len).rename(columns={'Age':'NumTotalVisits'})

# <codecell>

def make_impaired_call(val):
    if val < 10:
        return 'Impaired'
    elif val >= 10:
        return 'Not Impaired'
    else:
        return np.nan

hivd_data = pd.DataFrame({'HIVD score':redcap_data['TMHDS'],
                          'HIVD.I':redcap_data['TMHDS'].map(make_impaired_call)})
days_data = redcap_data

# <codecell>

rename_dict = {'Latest CD4 count (cells/uL)':'Latest CD4 count',
               'Latest viral load': 'Latest Viral load',
               'Latest CD8 count (cells/uL)':'Latest CD8 count',
               'Nadir CD4 count (cells/uL)':'Nadir CD4 count',
               'Nadir CD4 count (cells/uL)':'Nadir CD8 count',
               'Peak viral load (copies/mL)':'Peak viral load',
               }
outredcap = pd.concat([redcap_data.rename(columns=rename_dict), at_sample, to_sample, all_sample, num_data,
                       pd.DataFrame({'HAART':haart_data}), hivd_data], 
                      axis=1)

# <codecell>

outredcap[col_order[2:]]

# <codecell>

pd.rolling_mean?

# <codecell>


