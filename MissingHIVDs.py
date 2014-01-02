# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from __future__ import division
import os, os.path
from pandas import *
import numpy as np
import csv
from itertools import groupby
from collections import defaultdict
from copy import deepcopy
from types import TupleType
import matplotlib.pyplot as plt
from itertools import product, imap, islice
from patsy import dmatrices
from patsy.contrasts import Treatment
import statsmodels.api as sm
import shutil
import glob
import shlex
from subprocess import check_call, CalledProcessError
from tempfile import mkdtemp
from concurrent.futures import ProcessPoolExecutor
from types import TupleType
from rpy2.robjects import Formula
import rpy2.robjects as robjects
import pandas.rpy.common as com
from rpy2.rinterface import RRuntimeError
import sys


sys.path.append('/home/will/PySeqUtils/')
os.chdir('/home/will/HIVSystemsBio/NewCytokineAnalysis/')

import Rtools
import PlottingTools

# <codecell>

store_loc = os.path.join('home', 'will', 'HIVReportGen', 'Data',
                        'SplitRedcap', '2013-01-16',
                        'EntireCohort.hdf')
store = HDFStore('/'+store_loc)
redcap_data = store['redcap']

#fix the Visit Number issue
t = redcap_data['Event Name'].dropna().apply(lambda x: x.split(' - ')[0])
vnum_field = 'Patient visit number'
redcap_data[vnum_field] = redcap_data[vnum_field].combine_first(t)

# <codecell>

from dateutil.parser import parse
data = redcap_data[['Patient ID', vnum_field, 'Date of visit', 'Total Modified Hopkins Dementia Score']]
wrong_date = parse('2013-02-19 00:00:00')
data['Date of visit'][data['Date of visit'] == wrong_date] = np.nan

# <codecell>

first_score_date = min(data.dropna()['Date of visit'])
wanted_dates = (data['Date of visit'] > first_score_date)
check_visits = after_start & data['Total Modified Hopkins Dementia Score'].isnull()
tdata = data[check_visits | data['Date of visit'].isnull()]
tdata.to_excel('missing_hivd_data.xlsx', index = False)

# <codecell>

first_score_date

# <codecell>

data['Total Modified Hopkins Dementia Score'][wanted_dates | data['Date of visit'].isnull()].notnull().mean()

# <codecell>

data['Date of visit'].isnull().sum()

# <codecell>


