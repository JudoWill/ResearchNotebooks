# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from __future__ import division
import numpy as np
from functools import partial
from scipy.stats import linregress
import matplotlib.pyplot as plt
from matplotlib import dates
import pandas as pd
import gspread
from StringIO import StringIO
import csv
import sys
sys.path.append('/home/will/PatientPicker/')
sys.path.append('/home/will/PySeqUtils/')
import LoadingTools
from GeneralSeqTools import fasta_reader, fasta_writer, seq_align_to_ref

# <codecell>

pat_data = LoadingTools.load_redcap_data().groupby(['Patient ID', 'VisitNum']).first()

# <codecell>

mask = pat_data['Current Tobacco Use'] == 'No'
pat_data['Tobacco Use (packs/year)'][mask] = pat_data['Tobacco Use (packs/year)'][mask].fillna(0)
ewms = partial(pd.ewma, span = 2)
pat_data['Smoothed-Tobacco-Use'] = pat_data['Tobacco Use (packs/year)'].groupby(level = 'Patient ID').transform(ewms)

# <codecell>

cols = ['Tobacco Use (packs/year)', 'Smoothed-Tobacco-Use','Date Of Visit']
grouper = pat_data[cols].dropna().groupby(level = 'Patient ID')
fig, axs = plt.subplots(2, 1, figsize = (10, 10))
for pat, group in grouper: 
    if len(group) > 3:
        axs[0].plot_date(group['Date Of Visit'], group['Tobacco Use (packs/year)'], '-')
        axs[1].plot_date(group['Date Of Visit'], group['Smoothed-Tobacco-Use'], '-')
axs[0].set_ylabel('packs/year')
axs[1].set_ylabel('Smoothed packs/year')

# <codecell>

from sklearn.covariance import EllipticEnvelope
cytos = sorted(['IL.8','VEGF','IL.1beta',
        'G.CSF','EGF','IL.10','HGF',
        'FGF.basic','IFN.alpha','IL.6',
        'IL.12','Rantes','Eotaxin',
        'GM.CSF','MIP.1beta',
        'MCP.1','IL.5','IL.13', 'IFN.gamma','TNF.alpha',
        'IL.RA','IL.2','IL.7','IP.10',
        'IL.2R','MIG','IL.4','IL.15',
        'IL.17','MIP.1alpha']) + ['Th1', 'Th2']

cyto_data_raw = pd.read_csv('/home/will/HIVSystemsBio/NewCytokineAnalysis/CytoRawData.csv', sep = '\t')
cyto_data_raw['Th1'] = cyto_data_raw['IFN.gamma'] + \
                            cyto_data_raw['IL.2']+cyto_data_raw['TNF.alpha']
cyto_data_raw['Th2'] = cyto_data_raw['IL.4'] + \
                            cyto_data_raw['IL.5']+cyto_data_raw['IL.10']

# <codecell>

cyto_data = cyto_data_raw.groupby(['Patient ID', 'VisitNum']).mean()
tranfer_cols = ['Log-Latest-VL', 
                'Keep',
                'IsMale',
                'Race-Black',
                'Age',
                'HAART-Naive',
                'HAART-Non-Adherent',
                'HAART-Off',
                'HAART-On',
                'Hepatitis C status (HCV)']
for col in tranfer_cols:
    _, cyto_data[col] = cyto_data.align(pat_data[col], join='left', axis = 0)
cyto_data['HCV'] = cyto_data['Hepatitis C status (HCV)']

# <codecell>


