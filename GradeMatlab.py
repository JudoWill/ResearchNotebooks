# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import os, os.path
from subprocess import check_output
import shlex
import glob
from StringIO import StringIO
import contextlib
import tempfile
import shutil
import re


os.chdir('/home/will/matlabtmp/')

matlab_path = '/home/will/MATLAB/R2013a/bin/matlab'

# <codecell>

@contextlib.contextmanager
def tmp_changedir(new_dir):

    cur_dir = os.path.abspath(os.curdir)
    try:
        os.chdir(new_dir)
        yield new_dir
    finally:
        os.chdir(cur_dir)

# <codecell>

def get_warnings(mlint_out):
    
    last_line = 'For product information, visit www.mathworks.com.'
    handle = StringIO(mlint_out)
    for line in handle:
        if last_line in line:
            break
    return [line.strip() for line in handle if line.strip()]

def mlint_code(fpath):
    mlint_path = '/home/will/.matlab/R2013a/MLintGradingSettings.txt'
    mlint_cmd = "disp(checkcode('%(filename)s', '-config=%(gpath)s', '-notok', '-string')); quit;"
    test_mlint = mlint_cmd % {'filename':fpath, 'gpath':mlint_path }
    cmd = shlex.split(matlab_path + ' -nodesktop -nosplash -r "' + test_mlint + '"')
    out = check_output(cmd)
    return get_warnings(out)

# <codecell>

def grade_homework(grade_tmp, fpath, new_fname, prefix = None, req_files = [], sanitize_func = None):
    fname = fpath.rsplit(os.sep, 1)[-1]
    tmpdir = tempfile.mkdtemp(prefix=prefix, dir = '/home/will/matlabtmp/')
    grade_tmp_name = grade_tmp.split(os.sep)[-1].split('.')[0]
    with open(os.path.join(tmpdir, grade_tmp_name+'.m'), 'w') as handle:
        with open(grade_tmp) as ihandle:
            tdict = {'filename':new_fname.split('.')[0]}
            handle.write(ihandle.read() % tdict)
    
    
    shutil.copyfile(fpath, os.path.join(tmpdir, new_fname))
    if sanitize_func is not None:
        sanitize_func(os.path.join(tmpdir, new_fname))
    for f in req_files:
        t_fname = f.split(os.sep)[-1]
        shutil.copyfile(f, os.path.join(tmpdir, t_fname))
        
    with tmp_changedir(tmpdir):
        mlint_lines = mlint_code(new_fname)
        cmd = matlab_path + ' -nodesktop -nosplash -r "%s; quit;"' % grade_tmp_name
        cmd_list = shlex.split(cmd)
        out = check_output(cmd_list)
    try:
        outres = float(re.findall('THIS IS THE GRADE!!!:([\-\.\d]*)', out)[0])
        wrong_lines = re.findall('THIS QUESTION WAS WRONG!: .*', out)
        if (outres == 0) and (len(wrong_lines) == 0):
            wrong_lines = ["DIDN'T RUN!!"]
        shutil.rmtree(tmpdir)
        return mlint_lines, wrong_lines, outres
    except:
        return mlint_lines, ["DIDN'T RUN!!"], None

# <codecell>

def process_filenames(path):
    
    fname = path.rsplit(os.sep,1)[-1]
    parts = fname.split('_',4)
    bwid = parts[1]
    subname = parts[-1]
    date = parts[3].replace('_','-')
    
    return path, fname, date, bwid, subname

# <codecell>

def get_fnames(base_dir):
    
    for root, dirs, files in os.walk(base_dir):
        for f in sorted(files):
            if f.endswith('.m'):
                yield os.path.join(root, f)

# <codecell>

from functools import partial

def fix_micro_test_name(infile):
    with open(infile) as handle:
        tdata = handle.read()
    tdata = tdata.replace('microarray_test', 'MICROARRAY_TEST')
    with open(infile, 'w') as handle:
        handle.write(tdata)

dropbox_path = '/home/will/Dropbox/BMES375/Summer13/'
hw1_grade = partial(grade_homework, 
                    dropbox_path+'HW1/GradeHW1.m', 
                    req_files = [dropbox_path+'HW1/PAT_DATA.mat'])
hw2_grade = partial(grade_homework, 
                    dropbox_path+'HW2/GradeHW2.m', 
                    req_files = [dropbox_path+'HW2/CORRECT_DATA.mat',
                                 dropbox_path+'HW2/migraineSymptoms.tsv',
                                 dropbox_path+'HW2/migraineTreats.tsv',
                                 ])
micro_mid_grade = partial(grade_homework, 
                    dropbox_path+'MICRO_TEST/GradeMicroMid.m', 
                    req_files = [dropbox_path+'MICRO_TEST/MICRO_SOL.mat',
                                 ],
                    sanitize_func = fix_micro_test_name)
pred_mid_grade = partial(grade_homework, 
                    dropbox_path+'PRED_TEST/GradePredTest.m', 
                    req_files = [dropbox_path+'PRED_TEST/HW5_data.mat',
                                 dropbox_path+'PRED_TEST/COR_ANSWERS.mat',
                                 ],
                    sanitize_func = fix_micro_test_name)



grade_checks = [
                {
                 'name':'PredTest',
                 'grade_func':pred_mid_grade,
                 'grade_path':dropbox_path+'PRED_TEST/PRED_TEST_grades.csv',
                 'grade_files':partial(get_fnames, dropbox_path+'PRED_TEST/downloaded/'),
                 'out_path':dropbox_path+'PRED_TEST/PRED_TEST_graded/'
                 },
                #{
                # 'name':'MicroTest',
                # 'grade_func':micro_mid_grade,
                # 'grade_path':dropbox_path+'MICRO_TEST/MICRO_TEST_grades.csv',
                # 'grade_files':partial(get_fnames, dropbox_path+'MICRO_TEST/downloaded/'),
                # 'out_path':dropbox_path+'MICRO_TEST/MICRO_TEST_graded/'
                # },
                #{
                # 'name':'HW2',
                # 'grade_func':hw2_grade,
                # 'grade_path':dropbox_path+'HW2/HW2_grades.csv',
                # 'grade_files':partial(get_fnames, dropbox_path+'HW2/downloaded/'),
                # 'out_path':dropbox_path+'HW2/HW2_graded/'
                # },
                #{
                # 'name':'HW1',
                # 'grade_func':hw1_grade,
                # 'grade_path':dropbox_path+'HW1/HW1_grades.csv',
                # 'grade_files':partial(get_fnames, dropbox_path+'HW1/downloaded/'),
                # 'out_path':dropbox_path+'HW1/HW1_graded/'
                # },
                
                ]

# <codecell>

def rename_file(path, subname, count):
    if count > 0:
        tname = subname.split('.')[0] + '_' + str(count) + '.m'
    else:
        tname = subname.split('.')[0] + '.m'
    
    return os.path.join(path, tname)

# <codecell>

import csv
from itertools import islice

for row in grade_checks:
    grade_func = row['grade_func']
    with open(row['grade_path'], 'a') as handle:
        writer = csv.writer(handle)
        for f in row['grade_files']():
            f, fname, date, bwid, subname = process_filenames(f)
            mlint_lines, grade_lines, grade = grade_func(f, subname, prefix = bwid+'_')
            
            count = 0
            new_name = rename_file(row['out_path'], subname, count)
            while os.path.exists(new_name):
                count += 1
                new_name = rename_file(row['out_path'], subname, count)
            
            nf = new_name.split('/')[-1]
            print row['name'], bwid, nf, date, len(mlint_lines), grade
            print grade_lines
            tup = (bwid, date, nf, grade, len(mlint_lines), 
                             '; '.join(m.strip() for m in mlint_lines), 
                            '; '.join(m.strip() for m in grade_lines), 
                            )
            writer.writerow(tup)
            
            
            if (grade is not None) and (grade > 0):
                shutil.move(f, new_name)
            else:
                shutil.copy(f, new_name)
                

# <codecell>

.

