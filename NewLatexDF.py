# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from __future__ import division
from pandas import *
import numpy as np
import pickle
import os, os.path
from StringIO import StringIO
import glob
os.chdir('/home/will/HIVReportGen/tmp/')

# <codecell>

datasets = [pickle.load(open(f)) for f in glob.glob('example_pkls/*.pkl')]

# <codecell>

tmp_data = datasets[0]

# <codecell>

def process_index(index, other_names = None):
    nrows = len(index.names)
    if other_names:
        tmp = []
        for n, name in enumerate(other_names):
            if nrows == 1:
                index = index.insert(n, name)
            else:
                nones = (None,)*(nrows-1)
                nones += (name,)
                index = index.insert(n, nones)
        #print index
    
    if nrows == 1:
        return [list(index)], 1
    
    new_cols = []
    last_head = None
    for head_col in index:
        if last_head is None:
            new_cols.append(head_col)
            last_head = head_col
        else:
            new_head = []
            for t_item, l_item in zip(head_col, last_head):
                if t_item != l_item:
                    new_head.append(t_item)
                else:
                    new_head.append('')
            new_cols.append(new_head)
            last_head = head_col
    return map(list, zip(*new_cols)), nrows

# <codecell>

from types import StringType, IntType, FloatType

def treat_val(val):
    tp = type(val)
    if val is None:
        return ''
    elif (tp == StringType) or (tp == np.string_):
        #rint 'is string', val
        return val.strip()
    elif (tp == IntType) or (int(val) == val):
        return str(int(val))
    elif (tp == FloatType) or (tp == np.float64):
        if val < 1: #likely a p-value
            return '%.3E' % val
        else:
            return str(np.around(val, decimals=2))
    return str(val)


def process_dataset(df):
    
    drop_mask = (df==0).all(axis = 1)
    inds = drop_mask[drop_mask].index
    df = df.drop(inds, axis = 0)
    
    column_headers, depth = process_index(df.columns, other_names = df.index.names)
    
    row_index, ind_depth = process_index(df.index)
    
    width_per_char = 0.2 #cm/char
    index_width_max = 5 #cm
    row_index_widths = [max(len(x)*width_per_char for x in row) for row in row_index]
    row_capped_widths = [min()]
    
    out_buf = StringIO()
    ncols = len(column_headers[0])
    col_order = ''.join(['p{%.2f*\\textwidth}' % n for n in norm_size])+'l'*(ncols-ind_depth)
    out_buf.write('\\begin{tabular}{%s}\n' % col_order)
    out_buf.write('\\toprule\n')
    for row_num in range(depth):
        for v in column_headers[row_num][:-1]:
            out_buf.write(treat_val(v) + ' & ')
        out_buf.write(treat_val(column_headers[row_num][-1]) + '\\\\ \n')
    out_buf.write('\\midrule\n')
    
    tdata = df.values
    ind_data = np.array(row_index).transpose()
    
    nrows, ncols = tdata.shape
    
    for row in range(nrows):
        for col in range(ind_depth):
            out_buf.write(treat_val(ind_data[row, col]) + ' & ')
        
        for col in range(ncols-1):
            out_buf.write(treat_val(tdata[row, col]) + ' & ')
        out_buf.write(treat_val(tdata[row, ncols-1]) + ' \\\\ \n')
    
    
    
    out_buf.write('\\bottomrule\n')
    out_buf.write('\\end{tabular}\n')
    out_buf.seek(0)
    return out_buf.read()
    
    
    
with open('tmp.tex', 'w') as handle:
    for df in datasets:
        #print df
        out = process_dataset(df)
        handle.write(out)
        handle.write('\\newpage\n\n\n')

# <codecell>

from scipy.optimize import minimize
def fix_runs(in_val):
    
    mask = np.concatenate([in_val[:1] != in_val[:1], in_val[:-1] == in_val[1:]])
    out_val = in_val.copy()
    out_val[mask] = ''
    return out_val

def check_fun(widths, width_mat):
    txt_width = 12 #cm
    
    extra = widths - width_mat
    badness = -extra[extra<0].sum()
    if widths.sum() > txt_width:
        t = (widths.sum() - txt_width)*width_mat.shape[0]
        if t < 0:
            raise TypeError
        badness += (widths.sum() - txt_width)*width_mat.shape[0]
    
    return badness
    


def new_layout_table(df):
    
    
    cm_per_char = 0.5
    
    drop_mask = (df==0).all(axis = 1)
    inds = drop_mask[drop_mask].index
    df = df.drop(inds, axis = 0).dropna(axis = 1, how = 'all')
    
    header_rows = len(df.columns.names)
    index_cols = len(df.index.names)
    
    rep_df = df.reset_index().T.reset_index().T
    
    for row in range(header_rows):
        rep_df.values[row, :] = fix_runs(rep_df.values[row, :])
        
    for col in range(index_cols):
        rep_df.values[:, col] = fix_runs(rep_df.values[:, col])
    
    treated_df = rep_df.applymap(treat_val)
    widths = treated_df.applymap(len)
    req_cm_widths = widths*cm_per_char
    ncols = widths.values.shape[1]
    guess = (14/ncols)*np.ones((1,ncols))
    res = minimize(check_fun, guess, method = 'SLSQP', args = (req_cm_widths.values,), bounds = [(1, 12) for _ in range(ncols)])
    col_widths = list(res.x)
    
    out_buf = StringIO()
    
    col_order = ''.join(['p{%.2fcm}' % n for n in col_widths])
    out_buf.write('\\begin{tabular}{%s}\n' % col_order)
    out_buf.write('\\toprule\n')
    
    tdata = treated_df.values
    for row in range(header_rows):
        for col in range(ncols - 1):
            out_buf.write(tdata[row, col] + ' & ')
        out_buf.write(tdata[row, col+1] + ' \\\\ \n')
    out_buf.write('\\midrule \n')
    for row in range(header_rows, tdata.shape[0]):
        for col in range(ncols - 1):
            out_buf.write(tdata[row, col] + ' & ')
        out_buf.write(tdata[row, col+1] + ' \\\\ \n')
    
    
    
    out_buf.write('\\bottomrule\n')
    out_buf.write('\\end{tabular}\n')
    out_buf.seek(0)
    return out_buf.read()

with open('tmp.tex', 'w') as handle:
    handle.write("""\documentclass{article}
\usepackage{booktabs}
\\begin{document}
""")
    for df in datasets:
        w = new_layout_table(df)
        handle.write(w)
        handle.write('\\newpage')
    handle.write('\\end{document}')
    

# <codecell>

ndf = datasets[0]

# <codecell>

ndf.to_csv?

# <codecell>


