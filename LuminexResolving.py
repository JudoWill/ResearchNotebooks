# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import sys
sys.path.append('/home/will/PatientPicker/')

# <codecell>

import LoadingTools

# <codecell>

redcap_data = LoadingTools.load_redcap_data().groupby(['Patient ID', 'VisitNum']).first()

# <codecell>

cols = ['Date Of Visit']+[col for col in redcap_data.columns if col.startswith('Test-')]
cols

# <codecell>

import pandas as pd
pairs = [['A0138', 'A0360'],
         ['A0017', 'A0054'],
         ['A0235', 'A0266'],
         ['A0023', 'A0114'],
         ['A0059', 'A0403'],
         ]

res = []
for gp in pairs:
    tdata = redcap_data[cols].ix[gp].reset_index()
    tdata['TruePat'] = gp[0]
    res.append(tdata.copy())
    
nres = pd.concat(res, axis=0, ignore_index=True)

# <codecell>

nres.to_excel('/home/will/Downloads/luminex_resolve.xlsx')

# <codecell>


