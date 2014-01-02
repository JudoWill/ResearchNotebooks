# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import os, os.path
import re
from dateutil import parser
from pandas import *
import fileinput
import glob

os.chdir('/home/will/HIVTropism/')

# <codecell>

def yield_completion_lines():
    checker = re.compile('([\d:\-\s]*),\d*: INFO/MainProcess] Task ([\w\d.]*).*? succeeded in ([\d.]*)')
    handle = fileinput.FileInput(glob.glob('/home/will/celerylogs/*.log'))
    checks = (checker.findall(line) for line in handle)
    return [row[0] for row in checks if row]
   
def yield_task_lines():
    checker = re.compile('([\d:\-\s]*),\d*: INFO/MainProcess] Got task from broker: ([\w\d.]*).*?')
    handle = fileinput.FileInput(glob.glob('/home/will/celerylogs/*.log'))
    checks = (checker.findall(line) for line in handle)
    return [row[0] for row in checks if row]

# <codecell>

log_df = DataFrame(yield_completion_lines(), columns=['Time', 'Task', 'CompletionTime'])
log_df['CompletionTime'] = log_df['CompletionTime'].map(float)
log_df['Time'] = log_df['Time'].map(parser.parse)
grouped_data = log_df.groupby(['Time', 'Task'], as_index = False).mean()
pivoted_data = grouped_data.pivot(index = 'Time', columns='Task', values = 'CompletionTime')
completion_times = pivoted_data.resample('1Min', how = 'mean')

started_df = DataFrame(yield_task_lines(), columns = ['Time', 'Task'])
started_df['Time'] = started_df['Time'].map(parser.parse)
start_times = crosstab(started_df['Time'], started_df['Task']).resample('1Min', how = 'sum')

rolling_mean(completion_times, window = 30, min_periods = 1).plot(use_index = True, figsize=(10,3), logy=True, legend = False)
plt.ylabel('Completion Times (s)');
expanding_sum(start_times).plot(use_index = True, figsize=(10,3), legend = False)
plt.ylabel('Tasks Started');
#started_df.groupby(['Time', 'Task'], as_index = False).count()#.pivot(index = 'Time', columns = 'Task')

rolling_sum(start_times, window = 30, min_periods=1).plot(use_index = True, figsize=(10,3), legend = False)
plt.ylabel('Tasks Started per 30 min');

# <codecell>

log_df.tail()

# <codecell>


