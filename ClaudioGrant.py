# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from statsmodels.stats.power import tt_ind_solve_power

# <codecell>


normmean = 40
normstd = 17
admean = 55
adstd = 20

ef = abs(normmean-admean)/(adstd+normstd)
ratio = 1.0
alpha = 0.05
power = 0.9

n0 = tt_ind_solve_power(effect_size=ef, alpha = alpha, power = power, ratio = ratio)
print n0

# <codecell>

import pandas as pd

# <codecell>

data = pd.read_excel('/home/will/ClaudioStuff/TableData.xlsx', 'Sheet1')

# <codecell>

data['HIV'] = data['HIVInfected'] == 'pos'
data['Aged'] = data['Age']>50

# <codecell>

tab = pd.pivot_table(data,
                     cols = ['HIV', 'Impaired'],
                     rows = 'Aged',
                     values = 'ID',
                     aggfunc = 'count').fillna(0)

# <codecell>

tab.to_excel('/home/will/ClaudioStuff/count_table.xlsx')

# <codecell>

tab

# <codecell>


