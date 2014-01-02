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

os.chdir('/home/will/HIVSystemsBio/NewCytokineAnalysis/')

# <codecell>

cyto_data = read_csv('CytoRawData.csv', sep = '\t')

# <codecell>

import rpy2.robjects as robjects
import pandas.rpy.common as com

# <codecell>

r_dataframe = com.convert_to_r_dataframe(cyto_data)



robjects.r("""f <- function(r) {
            print(names(r))
        }      
        """)

# <codecell>

robjects.r['f'](r_dataframe)

# <codecell>

robjects.r('require(preprocessCore)')

# <codecell>

robjects.r("""quantnorm <- function(inputmatrix)
{
y<-normalize.quantiles(inputmatrix)
return(y)
}""" )

# <codecell>

small_dataframe = com.convert_to_r_matrix(cyto_data[['VEGF', 'HGF', 'Rantes']].T)

output = robjects.r['quantnorm'](small_dataframe)

# <codecell>

normed = com.convert_robj(output).T

# <codecell>

def quantile_norm_with_R(input_df):
    
    R_norm_func = robjects.r("""quantnorm <- function(inputmatrix)
{
y<-normalize.quantiles(inputmatrix)
return(y)
}""" )
    
    R_matrix = com.convert_to_r_matrix(input_df)
    print input_df
    normed_matrix = R_norm_func(R_matrix)
    normed_df = com.convert_robj(normed_matrix)
    print normed_df
    normed_df.index = input_df.index
    normed_df.columns = input_df.columns
    
    return normed_df
    
    
    
    

# <codecell>

norm_cyto_data = read_csv('CytoRawDataNorm.csv', sep = '\t')
known_pat_data = read_csv('CytoPatData.csv', index_col=[0,1], sep = '\t')

agg_cyto_data = norm_cyto_data.groupby(['Patient ID', 'VisitNum']).agg('median')
pat_cyto_data = merge(known_pat_data, agg_cyto_data,
                       left_index = True, right_index = True,
                        how = 'outer')

cytos = sorted(['IL.8','VEGF','IL.1beta',
        'G.CSF','EGF','IL.10','HGF',
        'FGF.basic','IFN.alpha','IL.6',
        'IL.12','Rantes','Eotaxin',
        'GM.CSF','MIP.1beta',
        'MCP.1','IL.5','IL.13', 'IFN.gamma','TNF.alpha',
        'IL.RA','IL.2','IL.7','IP.10',
        'IL.2R','MIG','IL.4','IL.15',
        'IL.17','MIP.1alpha'])

# <codecell>

robjects.r('require(nlme)')

# <codecell>

tmp = pat_cyto_data[['Gender', 'Age', 'Race', 'IL.12']]
tmp_R_df = com.convert_to_r_dataframe(tmp.reset_index())

# <codecell>

relevel_fun = robjects.r("""relevel_obj <- function(inputframe, column, ref)
{
print(inputframe[,column])

inputframe[,column] <- relevel(factor(inputframe[,column]), ref = ref)
return(inputframe)
}""")

# <codecell>

tmp_R_df = relevel_fun(tmp_R_df, 'Race', 'White')
tmp_R_df = relevel_fun(tmp_R_df, 'Gender', 'M')

# <codecell>

str(tmp_R_df.rx2('Age').nlevels)

# <codecell>

lme_model_func = robjects.r("""lmefunc <- function(eqn, data, reqn){
lm1 <- lme(eqn, random = reqn, data = data, na.action = na.omit)
return(summary(lm1))
}""")

# <codecell>

from rpy2.robjects import Formula

eqn = Formula('IL.12 ~ as.factor(Race) + as.factor(Gender) + Age')
rand_eqn = Formula('~1|Patient.ID')

summary = lme_model_func(eqn, tmp_R_df, rand_eqn)

# <codecell>

names = summary.names

# <codecell>

tdict = {}
for name in names:
    try:
        tdict[name] = com.convert_robj(summary.rx2(name))
    except:
        print name
    #except Exception:
    #    print name

# <codecell>

tdict['tTable']

# <codecell>

name

# <codecell>


