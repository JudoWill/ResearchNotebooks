# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <headingcell level=1>

# Are We On Drugs

# <markdowncell>

# I want to look at wether the patients Total Modified Hopkins Dementia Score is altered by whether the patient tested positive AT the visit. There has been some nominal concern that patients who are testing positive for Cocaine (or anything else) are still high (or recovering) and that could artificially lower thier TMHDS. The doctors in the clinic should be sending innebriated people home instead of testing them, but this is a nice second check.
# 
# So to look for this effect I am finding a set of patients in which they tested positive at one (or more visits) and then tested negative at one (or more visits) and had TMHDS at the corresponding visits. Using a Wilcoxen test (since the scores are non-normal) I'll look to see whether the visits with a positive drug test are lower then visits with a negative drug test.

# <codecell>

import os, os.path
from pandas import *
from matplotlib import pyplot as plt

# <codecell>

from pandas import HDFStore #This is a file storage format for large collections of data


store = HDFStore('/home/will/HIVReportGen/Data/BaseRedcap/HIVAIDSGeneticAnalys_DATA_LABELS_2013-01-16_1211.hdf')
redcap_data = store['redcap']
store.close()

# <rawcell>

# We need to pull out the columns with a drug-test, the TMHDS, and the patient visit information.

# <codecell>

drug_cols = [
'Amphetamines',
'Barbiturates',
'Benzodiazepines',
#'Cannabinoid', #remove Cannabus since the test can be positive for weeks after use.
'Cocaine + metabolite',
'Opiates',
'Phencyclidine'
]

data = redcap_data[drug_cols + ['Patient ID', 'Patient visit number', 'Total Modified Hopkins Dementia Score']]
data = data.rename(columns = {
                        'Patient visit number':'VisitNum',
                        'Total Modified Hopkins Dementia Score':'HIVDI'}).dropna()

# <rawcell>

# Look for rows with a positive drug test. Then group by Patient ID to calculate the mean of TMHDS with a positive test and the mean of the negative test. Then merge the datasets so we can find patients where we hav examples of both.

# <codecell>

pos_mask = data[drug_cols].any(axis = 1)
pos_scores = data[['Patient ID', 'VisitNum', 'HIVDI']][pos_mask]
neg_scores = data[['Patient ID', 'VisitNum', 'HIVDI']][~pos_mask]

pat_neg = neg_scores.groupby('Patient ID').mean()
pat_pos = pos_scores.groupby('Patient ID').mean()

merged_data = merge(pat_neg, pat_pos, 
                    left_index = True, 
                    right_index = True,
                    suffixes = ('_neg', '_pos'))

# <codecell>

merged_data.boxplot();
plt.ylabel('TMHDS');

# <rawcell>

# The boxplot shows the TMHDS for the SAME set of patients during a positive drug-test and a negitive drug test. 
# From looking at the Boxplot I don't see anything at all but I'll check anyways.

# <codecell>

from scipy.stats import wilcoxon

_, pval = wilcoxon(merged_data['HIVDI_neg'], merged_data['HIVDI_pos'])
print 'P-value:', pval
print 'Num Patients:', len(merged_data.index)

# <rawcell>

# Not even close! But this requires that a patient have BOTH a positive and a negative test. Our extreme drug users never have a negative test and are excluded from this test. I should also check to see whether the TMHDS is different between samples with a positive test vs. a negative test.

# <codecell>

data['Was-pos'] = pos_mask
data.boxplot(by = 'Was-pos');
plt.ylabel('TMHDS');

# <rawcell>

# In this boxplot I grouped all TMHDS from visits with a positive score (True-label) vs with a negative score (False label). 

# <codecell>

from scipy.stats import ks_2samp

_, pval = ks_2samp(pos_scores['HIVDI'], neg_scores['HIVDI'])
print 'P-value:', pval

# <markdowncell>

# Nope, I don't think the positive test alters the TMHDS. When I look at the same patient I don't see an effect of the positive test.

# <codecell>


