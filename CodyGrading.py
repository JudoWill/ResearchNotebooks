# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import csv
import urllib2
from urllib import urlencode
import os, os.path
from bs4 import BeautifulSoup
from dateutil.parser import parse
import re
import numpy as np
import sys
os.chdir('/home/will/Dropbox/BMES375/Summer13/')
sys.path.append('/home/will/PySeqUtils/')

# <codecell>

from PlottingTools import make_heatmap_df

# <codecell>

url_maps = {}
with open('cody_urls.csv') as handle:
    for row in csv.DictReader(handle):
        url = row['Cody Url']
        player_id = url.rsplit('/', 1)[-1].split('-')[0]
        url_maps[row['BannerWebID']] = player_id

# <codecell>

cody_problems = []
with open('cody_problems.csv') as handle:
    for row in csv.DictReader(handle):
        cody_problems.append((row['CodyProblem'], parse(row['DueDate'])))

# <codecell>


def get_web_data(user_num, problem_num):
    term = 'term='
    term += 'player_id%3A'+str(user_num)
    term += '+problem_id%3A'+str(problem_num)
    url = 'http://www.mathworks.com/matlabcentral/cody/solutions'
    return urllib2.urlopen(url + '?' + term).read()


def get_sols(soup):
    return soup.find_all('div', attrs = {'class':'solution-metric'})

def get_solution_size(solution):
    size = solution.span.text
    if size == 'Incorrect':
        return np.inf
    else:
        return int(size)
    
def get_grade(user, problem):
    soup = BeautifulSoup(get_web_data(user, problem))
    sols = get_sols(soup)
    for sol in sols:
        yield get_solution_size(sol)
    else:
        yield None

# <codecell>

from itertools import product
from concurrent.futures import ProcessPoolExecutor
from itertools import imap

def linker(tup):
    (bw_id, user_num), (problem, due_date) = tup
    tmp = []
    for size in get_grade(user_num, problem):
        tmp.append((bw_id, 'P-'+problem, size))
    return tmp

with ProcessPoolExecutor(max_workers = 10) as ex:
    tmp_data = []
    checks = product(url_maps.items(), cody_problems)
    for res in ex.map(linker, checks):
        tmp_data += res
        

# <codecell>

import pandas as pd

grade_df = pd.DataFrame(tmp_data, columns = ['BW_id', 'Prob', 'Size'])

# <codecell>

pdata = pd.pivot_table(grade_df, rows = ['BW_id', 'Prob'], 
                        values = 'Size', aggfunc = [np.min, len]).dropna()
pdata.head()

# <codecell>


best_scores  = pd.pivot_table(pdata.reset_index(), rows = 'BW_id', 
                                  cols = 'Prob', values = 'amin')
fig = make_heatmap_df(best_scores.T, figsize = (10,10))
plt.colorbar().set_label('Best Solution Nodesize')

# <codecell>

num_tries  = pd.pivot_table(pdata.reset_index(), rows = 'BW_id', 
                                  cols = 'Prob', values = 'len', 
                                  aggfunc = np.sum)
fig = make_heatmap_df(num_tries.applymap(np.log10).T, figsize = (10,10))
cb = plt.colorbar()
tickpos = [5, 10, 15, 20, 30, 40, 50, 60, 100]
cb.set_label('#Tries')
cb.set_ticks([np.log10(x) for x in tickpos])
cb.set_ticklabels(map(str, tickpos))

# <codecell>


