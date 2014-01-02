# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from Bio import Seq
from Bio import SeqIO
import pandas as pd
import numpy as np
import sys
import os
sys.path.append('/home/will/PySeqUtils/')
from GeneralSeqTools import fasta_reader, fasta_writer
from HIVAlignTools import SeqTransformer, build_aligners

# <codecell>

import HIVAlignTools

# <codecell>

HIVAlignTools.build_aligners()

# <codecell>

import shlex
from subprocess import check_call

def score_seq(known, guess, gapopen=10, gapextend=1):
    
    cmd = 'needle -asequence %(cb)s -bsequence %(seq)s -aformat score -gapopen %(go)f -gapextend %(ge)s -outfile %(out)s'
    with NamedTemporaryFile() as conb_handle:
        fasta_writer(conb_handle, [('SeqA', known)])
        conb_handle.flush()
        os.fsync(conb_handle.fileno())
        with NamedTemporaryFile() as seq_handle:
            fasta_writer(seq_handle, [('Seq1', guess)])
            seq_handle.flush()
            os.fsync(seq_handle.fileno())
            with NamedTemporaryFile() as out_handle:
                param_dict = {
                              'cb':conb_handle.name,
                              'seq':seq_handle.name,
                              'out':out_handle.name,
                              'go':gapopen,
                              'ge':gapextend
                              }
                cmd_list = shlex.split(cmd % param_dict)
                check_call(cmd_list)
                for line in out_handle:
                    parts = line.split()
                    if (len(parts) == 4):
                        return float(parts[-1][1:-2])
    


def score_seqs(known_seqs, guess_seqs, gapopen=10, gapextend=1):
    
    score = 0.0
    for ind in range(known_seqs.shape[0]):
        score += score_seq(known_seqs[ind], guess_seqs[ind],
                           gapopen=gapopen, gapextend=gapextend)
    return score



    

# <codecell>

from sklearn.base import BaseEstimator, ClusterMixin
from tempfile import NamedTemporaryFile
from Bio.SubsMat import MatrixInfo as matlist
from Bio.Blast import NCBIXML
from StringIO import StringIO
from Bio.Blast.Applications import NcbiblastxCommandline, NcbiblastnCommandline


class BlastAligner(BaseEstimator, ClusterMixin):
    
    def __init__(self, evalue=10, word_size=2, gapopen=11, gapextend=1, 
                 max_intron_length = 20, tmp_path = '/tmp/', result_type = 'aa',
                 db_path=NamedTemporaryFile(suffix='.fasta').name, num_threads=1):
        self.evalue = evalue
        self.word_size = word_size
        self.gapopen = gapopen
        self.gapextend = gapextend
        self.max_intron_length = max_intron_length
        self.tmp_path = tmp_path
        self.result_type = result_type
        self.db_path = db_path
        self.num_threads = num_threads
    
    def _write_seqs(self, X, handle):
        
        seqs = []
        for row in range(X.shape[0]):
            seq = ''.join(X[row])
            seqs.append(('Seq-%03i' % row, ''.join(l for l in seq if l.isalpha())))
            
        fasta_writer(handle, seqs)
        handle.flush()
        os.fsync(handle.fileno())
    
    
    def fit(self, X, y):
        
        
        empty_mask = y == 'XX'
        
        with open(self.db_path, 'w') as handle:
            self._write_seqs(y[~empty_mask], handle)
        cmd = 'makeblastdb -in %s -dbtype ' % self.db_path
        if self.result_type == 'aa':
            cmd += 'prot'
        else:
            cmd += 'nucl'
        
        check_call(shlex.split(cmd))
        
        
        return self
   
    
    def predict(self, X):
        
        if self.result_type == 'aa':
            blast_cmd = NcbiblastxCommandline
        else:
            blast_cmd = NcbiblastnCommandline
        
        
        with NamedTemporaryFile(dir=self.tmp_path, delete=True) as fasta_handle:    
            self._write_seqs(X, fasta_handle)
            blastx_cline = blast_cmd(query=fasta_handle.name,
                                     db = self.db_path, outfmt=5, 
                                     out = '-',
                                     evalue=self.evalue,
                                     word_size=self.word_size,
                                     gapopen=self.gapopen,
                                     gapextend=self.gapextend,
                                     max_intron_length=self.max_intron_length,
                                     num_threads=self.num_threads)
            stdout, stderr = blastx_cline()
        
        blast_records = NCBIXML.parse(StringIO(stdout))
        seqs = []
        names = []
        prots = []
        for rec in blast_records:
            for align in rec.alignments:
                hsp = align.hsps[0]
                prots.append({
                              'ID':rec.query,
                              'Seq':hsp.query
                              })
        blast_out = pd.DataFrame(prots).groupby('ID')['Seq'].first()
        wanted_out = pd.DataFrame({
                                   'ID':['Seq-%03i' % i for i in range(X.shape[0])],
                                   'want_seq':[True]*X.shape[0],
                                   }).groupby('ID')['want_seq'].first()
        out, _ = blast_out.align(wanted_out, join='right')
        
        return SeqTransformer().transform(out.fillna('XX').values)
    
    def score(self, X, y):
        
        empty_mask = y == 'XX'
        out_aligns = self.predict(X)
        
        pos_scores = score_seqs(y[~empty_mask], out_aligns[~empty_mask])
        bad_scores = score_seqs(out_aligns[empty_mask], out_aligns[empty_mask])
        return (pos_scores - bad_scores)/y.shape[0]
        

# <codecell>

neg_controls = {
                'env':['gag', 'pol', 'vif', 'vpr', 'ltr'],
                'gag':['ltr', 'vif', 'vpr', 'vpu', 'tat', 'rev', 'env'],
                #'ltr':['gag', 'pol', 'vpr', 'vpu', 'env'],
                'nef':['pol', 'gag', 'vpu', 'tat'],
                'pol':['env', 'vpr', 'vpu', 'nef', 'rev', 'ltr'],
                'rev':['ltr', 'gag', 'pol', 'vif', 'nef'],
                'tat':['ltr', 'pol', 'vif', 'nef'],
                'vif':['ltr', 'tat', 'vpu', 'rev', 'env', 'nef'],
                'vpr':['ltr', 'gag', 'pol', 'rev', 'env', 'nef'],
                'vpu':['ltr', 'gag', 'pol', 'vif', 'vpr', 'nef'],
                }

# <codecell>

def get_seq(prot_name, typ):
    trans_path = '/home/will/PySeqUtils/TransToolStuff/'
    tmp = 'HIV1_ALL_2012_%s_%s.fasta' % (prot_name.lower(), typ.upper())
    with open(trans_path + tmp) as handle:
        return SeqTransformer.get_from_fasta_handle(handle)


pos_names, pos_X = get_seq('genome', 'DNA')
env_names, env_y = get_seq('env', 'pro')
neg_names = []
neg_X = None
for neg_prot in neg_controls['env']:
    tnames, tx = get_seq(neg_prot, 'DNA')
    neg_names += tnames
    if neg_X is None:
        neg_X = tx.copy()
    else:
        neg_X = np.concatenate((neg_X, tx))
    
    

# <codecell>

pos_X_ser = pd.Series(pos_X, index=pos_names)
env_y_ser = pd.Series(env_y, index=env_names)

X_ser, y_ser = pos_X_ser.align(env_y_ser, join='inner')
X = X_ser.values
y = y_ser.values
in_env = set(env_names)
neg_inds = [num for num, name in enumerate(neg_names) if name not in in_env]
wneg_X = neg_X[neg_inds]
wneg_y = np.array(['XX']*wneg_X.shape[0])

print X.shape, y.shape, wneg_X.shape, wneg_y.shape

# <codecell>

Xall = np.concatenate((X, wneg_X))
yall = np.concatenate((y, wneg_y))

yclass = yall == 'XX'

# <codecell>

from sklearn.cross_validation import train_test_split, StratifiedKFold, cross_val_score, StratifiedShuffleSplit
from sklearn.grid_search import GridSearchCV

param_dict = {'evalue':np.logspace(-100, 1, 20)}

cv = StratifiedShuffleSplit(yclass, n_iter=3, test_size=500, train_size=100)
aligner = BlastAligner(num_threads=100)
gd = GridSearchCV(aligner, param_dict, refit=False, cv=cv, verbose=5)
gd.fit(Xall, yall)

# <codecell>

import 

# <codecell>


# <codecell>


