# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from pandas import *
import os, os.path
import matplotlib.pyplot as plt

os.chdir('/home/will/BloodBrainBarierProject/')

# <codecell>

data = read_csv('AdhesionDonerData_newer.csv', sep = ',')
print data.to_string()

# <codecell>

fig, axes = plt.subplots(1,2, figsize = (10,5), sharey = True)

for ax, ct in zip(axes.flatten(), ['CD3', 'CD14']):
    data[data['CellType'] == ct].boxplot(ax = ax)
    ax.set_title(ct)

# <codecell>

cor_vals

# <codecell>

cor_vals = read_csv('CorrData.csv', sep = '\t')
cor_vals['NGrouping'] = cor_vals['Grouping'].map(lambda x: x.split('-')[0])
cor_vals['Donor'] = cor_vals['Grouping'].map(lambda x: x.split('-')[1])
print cor_vals
print cor_vals['Grouping'].unique()

# <codecell>

tmp = cor_vals.groupby(['Plot Title', 'Grouping']).agg({'MFI':'mean', 'Number of adhering cells':'mean', 'NGrouping':'first', 'Donor':'first'})
print tmp.head()

# <codecell>

clist = dict([('untreated', 'o'),
         ('m24', '+'),
         ('m48', '*'),
         ('m72', 'D'),
        ('1', 'r'),
        ('2', 'g'),
        ('3', 'b')])

plots = sorted(cor_vals['Plot Title'].unique())
results = []
for p in plots:
    plt.figure(figsize = (10,10))
    
    plt.hold(True)
    tdata = tmp.ix[p]
    for _, row in tdata.iterrows():
        
        plt.scatter(row['MFI'], row['Number of adhering cells'], 
                    marker = clist[row['NGrouping']], color = clist[row['Donor']],
                    s = 100)
    for donor, df in tdata.groupby('Donor'):
        m, b, rqs, pval, _ = linregress(df['MFI'], df['Number of adhering cells'])
        results.append((p, donor, m, b, rqs, pval))
        xpos = np.linspace(df['MFI'].min()*0.9, df['MFI'].max()*1.1, 10)
        ypos = m*xpos+b
        plt.plot(xpos, ypos, color = clist[donor])
    
    plt.title(p)
    fname = 'TrendLines-'  + p.replace(' ', '-') + '.png'
    plt.hold(False)
    
    plt.savefig(fname)

# <codecell>

from scipy.stats import linregress

res = linregress(tmp['Expressors'].values, tmp['Adhesor'].values)
print res

# <codecell>

resdf = DataFrame(results, columns = ['Anal', 'Donor', 'm', 'b', 'RSquared', 'Pval'])

# <codecell>

resdf.to_excel('trend_results.xls')

# <codecell>

four_point = cor_vals.groupby(['Plot Title', 'NGrouping']).mean()
plots = sorted(cor_vals['Plot Title'].unique())
four_res = []
for p in plots:
    m, b, rqs, pval, _ = linregress(four_point.ix[p]['MFI'], four_point.ix[p]['Number of adhering cells'])
    four_res.append((p, rqs, pval))
four_df = DataFrame(four_res, columns = ['Anal', 'R^2', 'Pval'])
print four_df
four_df.to_excel('FourPoints.xls')

# <codecell>


