# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from __future__ import division
import pandas as pd
import csv
import sys
import os, os.path
import numpy as np
sys.path.append('/home/will/PySeqUtils/')
import GeneralSeqTools

# <codecell>

print 2+2

# <codecell>

import glob
from Bio import SeqIO
from concurrent.futures import ProcessPoolExecutor
from itertools import islice, imap
import os, os.path
import csv

def get_gi_acc(fname):
    gb = fname.split('/')[-1].split('.')[0]
    with open(fname) as handle:
        for line in handle:
            if line.startswith('ACCESSION'):
                acc = line.strip().split()[-1]
                return gb, acc
    raise AssertionError



gi_to_acc_dict = {}

fname = '/home/will/WLAHDB_data/gi_to_acc.csv'
if os.path.exists(fname):
    with open(fname) as handle:
        for row in csv.reader(handle):
            gi_to_acc_dict[row[0]] = row[1]
else:
    gb_files = glob.glob('/home/will/WLAHDB_data/GenbankDL/*.gb')
    with open(fname, 'w') as handle:
        writer = csv.writer(handle)
        for num, (gbm, acc) in enumerate(imap(get_gi_acc, gb_files)):
            if (num == 100) or (num % 50000 == 0):
                print num
            gi_to_acc_dict[gbm] = acc
            writer.writerow((gbm, acc))


# <codecell>



files = [('C', sorted(glob.glob('/home/will/WLAHDB_data/SeqDump/C_*'))),
         ('B', sorted(glob.glob('/home/will/WLAHDB_data/SeqDump/B_*')))]
seqs = []
for sub, sfiles in files:
    for f in sfiles:
        with open(f) as handle:
            base_name = f.rsplit(os.sep,1)[1].rsplit('.',1)[0]
            prot = base_name.split('_')[1]
            for name, seq in GeneralSeqTools.fasta_reader(handle):
                seqs.append({
                             'Seq':seq,
                             'ID':gi_to_acc_dict[name],
                             'Prot':prot,
                             'Subtype':sub
                             })
            
seqdf = pd.DataFrame(seqs)

# <codecell>

pseqdf = pd.pivot_table(seqdf, 
                        rows = ['Subtype', 'ID'], 
                        cols = 'Prot', 
                        values = 'Seq', 
                        aggfunc = 'first')

# <codecell>

checker = pd.read_csv('/home/will/Downloads/hivseqdb/HIVSeqDB-MetaData.tsv', sep = '\t')
check_res = pd.merge(checker, pseqdf.ix['B'],
                     left_on = 'Accession',
                     right_index=True,
                     how = 'inner')
check_res
               

#pseqdf

# <codecell>

mat_dict = {'B':'x4r5', 'C':'subC'}

v3_res = []
for sub, df in pseqdf.groupby(level = 'Subtype'):
    v3_seqs = [(gi, seq) for (_, gi), seq in df['v3'].dropna().to_dict().items()]
    
    for num, row in enumerate(GeneralSeqTools.WebPSSM_V3_series(v3_seqs, matrix=mat_dict[sub])):
        if (num == 0) or (num == 100) or (num == 10) or (num % 5000 == 0):
            print num
        v3_res.append(row)
        v3_res[-1]['Subtype'] = sub
        

# <codecell>

def safe_mean(vals):
    if len(vals.dropna()):
        return np.mean(vals.map(float).dropna())
    else:
        return np.nan
v3_df = pd.DataFrame(v3_res).convert_objects()
x4_pred = v3_df.groupby(['Subtype', 'name'])['pred'].agg(safe_mean) > 0.5

# <codecell>

conb_ltr = """tGGAaGGGcTaaTttacTCcCaaaaaAgacAagAtATcCTTGAtcTGTGG
gtctAccAcAcaCAaGGCTacTTCCCtGAtTgGCAgAAcTACACAccAGG
gccaGGgatcAGaTatCCacTgaccTTTGGaTGgTGcTtcAAgcTAGTAC
CAgTtgAgCcAGagaaggTagAagagGccAatgaaggagagaacaacagc
tTGtTaCAcCCtatgagCCtgCATGGgatgGAgGAcccgGAgaaagAAGt
gtTagtgTGGAagTttGACAgccgccTaGcatttcatCAcatggCccgaG
AgctgcATCCggAgTactacaaggActGcTGACatcgagctttctacaaG
GGACTTTCCgCtgGGGACTTTccagggagGcGtggcctGGgcgggaCtgg
GGAgtggCgagCCCtcAGAtgcTgCATATAAGCAGCtGCttTttGccTGT
ACtGGgTCTCTCTggttaGaCCAGATCtGAGCctGGGAGcTCtctggcta
actagggaacccactgcttaagcctcaataaagcttgccttgagtgcttc
aagtagtgtgtgcccgtctgttgtgtgactctggtaactagagatccctc
agacccttttagtcagtgtggaaaatctctagca""".replace('\n', '').upper()

conc_ltr = """tggaagggttaatttactctaagaaaaggcaagagatccttgatttgtgg
gtctatcacacaCaAGgctactTcCCtGAtTGGCAaaacTacACACCgGG
aCCaGGggtcAgatacCCacTgACctttGGaTGGtgcTtcAAgctaGTaC
CAgTtgacCCaaggGaagtaGAagAggccaacgaaggaGAaaacaaCtGt
tTgcTaCAcCCtatgagcCagCatGGaatgGAggAtgaacacagagAAgt
atTaaagTGGaagtttgacagtcaccTaGcacgcagacacatggcCcGcg
aGctacATCcGGAgTatTAcAaagACTGCTGacACagaaGGgACTtTccg
ctgggactttccactgGGGcgttccaggaggtgtggtctGGGCGGgActG
GGAGTGGtcaaCCCtCaGatGctgCATATAAGCagcTGCTTTtcgcctgt
actgggtctctctaggtagaccagatctgagcctgggagctctctggcta
tctagggaacccactgcttaagcctcaataaagcttgccttgagtgctct
aagtagtgtgtgccctctgttttgactctggtaactagagatccctcaga
cccttttggtagtgaggaaatctctagca""".replace('\n', '').upper()


ltrs = {'B':conb_ltr, 'C':conc_ltr}

# <codecell>

ltr_aligns = []
for align_sub, ltr_seq in ltrs.items():
    for sub, df in pseqdf.groupby(level = 'Subtype'):
        ltr_seqs = [(gi, seq) for (_, gi), seq in df['ltr'].dropna().to_dict().items()]
        for num, (gi, align) in enumerate(GeneralSeqTools.seq_align_to_ref(ltr_seqs, ltr_seq)):
            if (num == 0) or (num == 100) or (num == 10) or (num % 1000 == 0):
                print align_sub, sub, num
            
            ltr_aligns.append({
                               'Subtype':sub,
                               'ID':gi,
                               'Align-Con'+align_sub:align
                               })

# <codecell>

ltr_align_df = pd.pivot_table(pd.DataFrame(ltr_aligns),
                              rows = ['Subtype', 'ID'],
                              values = ['Align-ConC', 'Align-ConB'],
                              aggfunc = 'first')

# <codecell>

tmp_seqs = []
for col in ltr_align_df.columns:
    res = ltr_align_df[col].apply(lambda x: pd.Series(list(x)))
    res.columns = pd.MultiIndex.from_tuples([(col, pos+1) for pos in res.columns],
                                          names = ['AlignType', 'OCol'])
    tmp_seqs.append(res.copy())
    
ltr_cols = pd.concat(tmp_seqs, axis = 1).replace('-', np.nan)
ltr_cols

# <codecell>

from scipy.stats import fisher_exact
from sklearn.metrics import adjusted_mutual_info_score
from functools import partial

def fishers_apply(group_series, seq_series):
    
    vgroup, vseq = group_series.dropna().align(seq_series.dropna(), join = 'inner')
    
    let = vseq.value_counts().idxmax()
    ftable = [[(vseq[vgroup]==let).sum(), (vseq[vgroup]!=let).sum()],
              [(vseq[~vgroup]==let).sum(), (vseq[~vgroup]!=let).sum()]]
    ratio, pval = fisher_exact(ftable)
    mi = max(adjusted_mutual_info_score(vgroup.values, vseq.values), 0)
    true_con = (vseq[vgroup]==let).mean()
    false_con = (vseq[~vgroup]==let).mean()
    
    ser = pd.Series([ratio, pval, mi, true_con, false_con], 
                    index = ['Ratio', 'Pval', 'NMI', 'TrueCons', 'FalseCons'])
    
    return ser

v3_fishers = partial(fishers_apply, x4_pred.ix['C'])
col_pvals = ltr_cols['Align-ConC'].ix['C'].apply(v3_fishers, axis=0).T

# <codecell>

from statsmodels.sandbox.stats.multicomp import multipletests
from scipy.stats import gaussian_kde

reject, adjpval, _, _ = multipletests(col_pvals.dropna()['Pval'], method = 'fdr_bh', alpha = 0.01)
pval_threshold = adjpval[reject].max()

# <codecell>

fig, axs = plt.subplots(4,1, figsize = (10,8), sharex=True)

col_pvals['LogP'] = -np.log10(col_pvals['Pval'])
col_pvals['NMI-P'] = 100*col_pvals['NMI']
plot_cols = [('LogP', '-log10(p-value)'),
             ('NMI-P', '% Explained Variation'),
             ('TrueCons', 'X4 Cons'),
             ('FalseCons', 'R5 Cons')]

for (col, label), ax in zip(plot_cols, axs.flatten()):

    ax.plot(col_pvals.index, col_pvals[col])
    ax.set_ylabel(label)
    
    if col == 'LogP':
        ax.hlines(-np.log10(pval_threshold), 0, max(col_pvals.index), 'r')
    elif 'Cons' in col:
        ax.set_ylim([0,1])
    
    if ax.is_first_row():
        ax.set_title('SubC X4/R5')
    elif ax.is_last_row():
        ax.set_xlabel('ConC')
ax.set_xlim([0, max(col_pvals.index)])
fig.tight_layout()
plt.savefig('/home/will/Dropbox/Wigdahl HIV Lab/SubCAnalysis/subC_x4_r5.png')
    


# <codecell>

ld = []
for f in glob.glob('/home/will/SubCData/LANLRes/*.txt'):
    ld.append(pd.read_csv(f, sep='\t', index_col=0))
    
lanl_data = pd.concat(ld, axis = 0, ignore_index=True)

# <codecell>

pd.crosstab(lanl_data['Georegion'], lanl_data['Subtype'])

# <codecell>

geo_regions = lanl_data.groupby('Accession')['Georegion'].first()

# <codecell>

def mni_apply(group_series, seq_series):
    
    vgroup, vseq = group_series.dropna().align(seq_series.dropna(), join = 'inner')
    
    let = vseq.value_counts().idxmax()
    cons = (vseq == let).mean()
    mi = max(adjusted_mutual_info_score(vgroup.values, vseq.values), 0)
    
    ser = pd.Series([mi, cons], 
                    index = ['NMI', 'Cons'])
    
    return ser

geo_mni = partial(mni_apply, geo_regions)
col_mni = ltr_cols['Align-ConC'].ix['C'].apply(geo_mni, axis=0).T



# <codecell>

fig, axs = plt.subplots(2,1, figsize = (10,8), sharex=True)


col_mni['NMI-P'] = 100*col_mni['NMI']
col_mni['Cons-P'] = 100*col_mni['Cons']
plot_cols = [('NMI-P', '% Explained Variation'),
             ('Cons-P', '% Conservation')]

for (col, label), ax in zip(plot_cols, axs.flatten()):

    ax.plot(col_mni.index, col_mni[col])
    ax.set_ylabel(label)
    ax.set_ylim([0,100])
    
    if ax.is_first_row():
        ax.set_title('SubC Geo')
    elif ax.is_last_row():
        ax.set_xlabel('ConC')
ax.set_xlim([0, max(col_mni.index)])
fig.tight_layout()
plt.savefig('/home/will/Dropbox/Wigdahl HIV Lab/SubCAnalysis/subC_geo.png')

# <codecell>

from Bio.Seq import Seq
from Bio import Motif
from StringIO import StringIO
from itertools import groupby
from operator import methodcaller

def yield_motifs():
    with open('/home/will/LTRtfAnalysis/Jaspar_PWMs.txt') as handle:
        for key, lines in groupby(handle, methodcaller('startswith', '>')):
            if key:
                name = lines.next().strip().split()[-1].lower()
            else:
                tmp = ''.join(lines)
                mot = Motif.read(StringIO(tmp), 'jaspar-pfm')
                yield name, mot
                yield name+'-R', mot.reverse_complement()
pwm_dict = {}
for num, (name, mot) in enumerate(yield_motifs()):
    if num % 100 == 0:
        print num
    pwm_dict[name] = mot
    
    
tmp = u"""
A 0  0 6 1 0 0 0 4 2 2 0 0 3 
C  1 1 1 0 5 6 4 1 0 0 0 3 5 5 4 0  
G  0 6 0 1 1 0 0 0 0 7 1 1 0 0 1 0 
T  6 0 0 0 1 1 3 5 7 0 0 0 0 2 2 4"""

pwm_dict['coup2'] = Motif.read(StringIO(tmp), 'jaspar-pfm')
pwm_dict['coup2-R'] = Motif.read(StringIO(tmp), 'jaspar-pfm').reverse_complement()

# <codecell>

'Bio.Motif._Motif.Motif' in str(type(pwm_dict['coup2']))


# <codecell>

from Bio.Alphabet import IUPAC

def score_seq(mot, seq):
    bseq = Seq(seq, alphabet=IUPAC.unambiguous_dna)
    scores = mot.scanPWM(bseq)
    for pos, score in enumerate(scores.flatten(),1):
        if ~np.isnan(score):
            tseq = seq[pos:pos+len(mot)]
            yield pos, tseq, score

def get_region(seq, reference, regions = None):
    
    if regions == None:
        regions = [(300, 400)]
    
    tmp_seqs = [('conc', reference),
                ('guess', seq)]
    aligned = dict(GeneralSeqTools.call_muscle(tmp_seqs))
    out = []
    for _, start, stop in regions:
        conc_pos = 0
        align_start = None
        for align_pos, l in enumerate(aligned['conc']):
            if l != '-':
                conc_pos += 1
            if conc_pos == start:
                align_start = align_pos
            if conc_pos == stop:
                align_stop = align_pos
                break
        yield seq[align_start:align_stop].replace('-', '')

            
def count_mots(mot, threshes, seq):
    
    check_seq = get_region(seq).next()
    
    if len(check_seq) > 0:
        bseq = Seq(check_seq, alphabet=IUPAC.unambiguous_dna)
        scores = mot.scanPWM(bseq)
        out = []
        for name, thresh in threshes:
            out.append((name, (scores>=thresh).sum()))
        out.append(('300to400Size', len(check_seq)))
    else:
        
        out = [(n, np.nan) for n, _ in thresh]
        out += [('300to400Size', np.nan)]
    return pd.Series(dict(out))

def check_max_scores(mots, regions, reference, seq):
    
    names = []
    out_seqs = []
    out_scores = []
    check_seqs = get_region(seq, reference, regions = regions)
    for (name, rstart, rstop), check_seq, mot in zip(regions, check_seqs, mots):
        names.append(name)
        if len(check_seq) < len(mot):
            out_scores.append(np.nan)
            out_seqs.append(np.nan)
        else:
            bseq = Seq(check_seq, alphabet=IUPAC.unambiguous_dna)
            scores = np.maximum(mot.scanPWM(bseq),
                                mot.reverse_complement().scanPWM(bseq))
            start_pos = np.argmax(scores)
        
            out_scores.append(np.max(scores))
            out_seqs.append(check_seq[start_pos:start_pos+len(mot)])
        
    
    inds = list(product(['Scores', 'Seqs'],
                        names))
    
    mi = pd.MultiIndex.from_tuples(inds)
    ser = pd.Series(out_scores + out_seqs, 
                    index = mi)
    
    return ser
    
    

# <codecell>

conC_regions = [('NFAT-III', 163-5, 181+5),
                ('NFAT-II', 182-5, 200+5),
                ('NFAT-I', 240-5, 250+5),
                ('CEBP', 198-5, 210+5),
                ('NFKB-IV', 332-10, 332+10),
                ('NFKB-III', 340-5, 349+5),
                ('NFKB-II', 353-5, 362+5),
                ('NFKB-I', 366-5, 375+5),
                ('SpIII', 379-5, 387+5),
                ('SPII', 389-5, 398+5),
                ('SPI', 400-5, 408+5)]
conC_mots = [pwm_dict['nfatc2'],
             pwm_dict['nfatc2'],
             pwm_dict['nfatc2'],
             pwm_dict['cebpa'],
             pwm_dict['nf-kappab'],
             pwm_dict['nf-kappab'],
             pwm_dict['nf-kappab'],
             pwm_dict['nf-kappab'],
             pwm_dict['sp1'],
             pwm_dict['sp1'],
             pwm_dict['sp1']
             ]
conC_tf_checks = partial(check_max_scores, 
                         conC_mots,
                         conC_regions,
                         conc_seqs['LTR'])

# <codecell>


subC_scores = {}
for key, seq in pseqdf.ix['C']['ltr'].dropna().to_dict().items():
    subC_scores[key] = conC_tf_checks(seq)
subC_scores = pd.DataFrame(subC_scores).T

# <codecell>

conB_regions = [('AP1-IV', 104-5, 104+15),
                ('AP1-III', 119-5, 119+15),
                ('AP1-II', 154-5, 154+15),
                ('AP1-I', 213-5, 213+15),
                ('CEBP-II', 280-5, 280+15),
                ('ETS', 304-5, 304+15),
                ('ATF/Creb', 329-5, 329+15),
                ('CEBP-I', 329-5, 329+15),
                ('NFKB-II', 349-5, 349+15),
                ('NFKB-I', 362-5, 362+15),
                ('SpIII', 376-5, 376+15),
                ('SPII', 387-5, 387+15),
                ('SPI', 398-5, 398+15),]
conB_mots = [pwm_dict['ap1'],
             pwm_dict['ap1'],
             pwm_dict['ap1'],
             pwm_dict['ap1'],
             pwm_dict['cebpa'],
             pwm_dict['ets1'],
             pwm_dict['creb1'],
             pwm_dict['cebpa'],
             pwm_dict['nf-kappab'],
             pwm_dict['nf-kappab'],
             pwm_dict['sp1'],
             pwm_dict['sp1'],
             pwm_dict['sp1']]
conB_tf_checks = partial(check_max_scores, 
                         conB_mots,
                         conB_regions,
                         conb_seqs['LTR'])

# <codecell>

subB_scores = {}
for key, seq in pseqdf.ix['B']['ltr'].dropna().to_dict().items():
    subB_scores[key] = conB_tf_checks(seq)
subB_scores = pd.DataFrame(subB_scores).T

# <codecell>

ltr_scores = pd.concat([subC_scores, subB_scores], axis=0)
ltr_scores

# <codecell>

weblogolib.LogoData.from_counts?

# <codecell>

import weblogolib
from weblogolib.colorscheme import nucleotide
from corebio.seq import unambiguous_dna_alphabet

def make_logo(fasta_path, png_path, title, start_pos, counts = None):
    
    if counts:
        mat_counts = np.array([counts[l] for l in 'ACGT']).transpose()
        data = weblogolib.LogoData.from_counts(unambiguous_dna_alphabet, mat_counts)
    else:
        with open(fasta_path) as handle:
            seqs = weblogolib.read_seq_data(handle)
            data = weblogolib.LogoData.from_seqs(seqs)
    
    options = weblogolib.LogoOptions()
    options.logo_title = title
    options.resolution = 500
    options.first_index = start_pos
    options.number_interval = 1
    options.rotate_numbers = True
    options.color_scheme = nucleotide
    fmt = weblogolib.LogoFormat(data, options)
    with open(png_path, 'w') as handle:
        weblogolib.png_formatter(data, fmt, handle)
    

# <codecell>

subtypes = ['C', 'B']
cols = ltr_scores['Seqs'].columns
results = []
for sub, col in product(subtypes, cols):
    
    fix_name = col.replace('/', '-')
    start_pos = 0
    if sub == 'B':
        for (name, start, _), mot in zip(conB_regions, conC_mots):
            if col == name:
                start_pos = start+5
                break
    else:
        for (name, start, _), mot in zip(conC_regions, conB_mots):
            if col == name:
                start_pos = start+5
                break
    
    pred_counts, mask = ltr_scores.align(x4_pred.ix[sub].dropna(), 
                                         join = 'inner', 
                                         axis = 0)
    pwm_name = '/home/will/SubCData/TFfasta/PWM-%s.png' % fix_name
    pwm_name_r = '/home/will/SubCData/TFfasta/PWM-%s-R.png' % fix_name
    make_logo(None,
              pwm_name,
              fix_name,
              1, counts = mot.counts)
    make_logo(None,
              pwm_name_r,
              fix_name,
              1, counts = mot.reverse_complement().counts)
    
    if mask.mean() > 0.5:
        r5mask, x4mask = (mask, ~mask)
    else:
        r5mask, x4mask = (~mask, mask)
    x4_name = '/home/will/SubCData/TFfasta/X4-%s-%s' % (fix_name, sub)
    r5_name = '/home/will/SubCData/TFfasta/R5-%s-%s' % (fix_name, sub)
    with open(x4_name+'.fasta', 'w') as handle:
        x4_seqs = pred_counts['Seqs'][col][x4mask].dropna().to_dict().items()
        x4_scores = pred_counts['Scores'][col][x4mask].astype(float).dropna()
        GeneralSeqTools.fasta_writer(handle, x4_seqs)
    if len(x4_seqs) == 0:
        os.remove(x4_name+'.fasta')
        continue
    make_logo(x4_name+'.fasta',
              x4_name+'.png',
              'X4-%s-%s' % (sub, col),
              start_pos)
    
    with open(r5_name+'.fasta', 'w') as handle:
        r5_seqs = pred_counts['Seqs'][col][r5mask].dropna().to_dict().items()
        r5_scores = pred_counts['Scores'][col][r5mask].astype(float).dropna()
        GeneralSeqTools.fasta_writer(handle, r5_seqs)
    make_logo(r5_name+'.fasta',
              r5_name+'.png',
              'R5-%s-%s' % (sub, col),
              start_pos)
    #print pred_sp_counts['Scores'][name][mask]
    r, pval = ttest_ind(x4_scores, r5_scores)
    results.append({
                    'Subtype': sub,
                    'X4Binding': x4_scores.mean(),
                    'R5Binding': r5_scores.mean(),
                    'NX4': len(x4_seqs),
                    'NR5': len(r5_seqs),
                    'pvalue': pval,
                    'TF': col
                    })
    

# <codecell>

pd.DataFrame(results).to_excel('/home/will/Dropbox/Wigdahl HIV Lab/SubCAnalysis/TF_binding_scores.xlsx')

# <codecell>

print pred_sp_counts['Seqs'][name][mask].dropna().head()
print pred_sp_counts['Seqs'][name][~mask].dropna().head()

# <codecell>

threshes = [('FPR-0.01', Motif.Thresholds.ScoreDistribution(nfkb_mot, precision = 50).threshold_fpr(0.01)),
            ('FPR-0.001', Motif.Thresholds.ScoreDistribution(nfkb_mot, precision = 50).threshold_fpr(0.001)),
            ('FPR-0.002', Motif.Thresholds.ScoreDistribution(nfkb_mot, precision = 50).threshold_fpr(0.002)),
            ('FPR-0.005', Motif.Thresholds.ScoreDistribution(nfkb_mot, precision = 50).threshold_fpr(0.005)),
            ]

# <codecell>

nfkb_count = partial(count_mots, nfkb_mot, threshes)
nf_counts = pseqdf.ix['C']['ltr'].dropna().apply(nfkb_count)
nf_counts['InsertSize'] = nf_counts['300to400Size']-100
#print nf_counts.value_counts()

# <codecell>

pred_nf_counts, pred_x4_nf = nf_counts.align(subc_pred_x4.dropna(), join = 'inner', axis = 0)
pred_nf_counts['IsX4']  = pred_x4_nf

# <codecell>

t.plot?

# <codecell>

fig, ax = plt.subplots(1,1, figsize = (10,5))
edges = np.arange(-10, 30, 5)
x4counts, _ = np.histogram(pred_nf_counts['InsertSize'][mask], bins = edges)
r5counts, _ = np.histogram(pred_nf_counts['InsertSize'][~mask], bins = edges)

t = pd.DataFrame({
                  'X4':x4counts/x4counts.sum(),
                  'R5':r5counts/r5counts.sum()
                  }, index = edges[1:])
(100*t).plot(kind = 'bar', ax=ax, color = 'rb', grid = False, legend = False)
ax.set_ylim([0, 100])
print x4counts.sum(), r5counts.sum()
ax.set_ylabel('% Sequences')
ax.set_xlabel('Indel Size')
fig.savefig('/home/will/Dropbox/Wigdahl HIV Lab/SubCAnalysis/insert_size.png', dpi = 300)

# <codecell>

t.plot?

# <codecell>

fig, ax = plt.subplots(1,1, figsize = (10,5))
edges = sorted(pred_nf_counts['InsertSize'].dropna().unique())
x4counts, _ = np.histogram(pred_nf_counts['InsertSize'][mask], bins = edges)
r5counts, _ = np.histogram(pred_nf_counts['InsertSize'][~mask], bins = edges)

t = pd.DataFrame({
                  'X4':x4counts/x4counts.sum(),
                  'R5':r5counts/r5counts.sum()
                  }, index = edges[1:])
(100*t).plot(ax=ax, color = 'rb', grid = False, legend = False, linewidth = 5, alpha = 0.8)
ax.set_ylim([0, 50])
print x4counts.sum(), r5counts.sum()
ax.set_ylabel('% Sequences')
ax.set_xlabel('Indel Size')
fig.savefig('/home/will/Dropbox/Wigdahl HIV Lab/SubCAnalysis/insert_size_line.png', dpi = 300)

# <codecell>

tcols = ['InsertSize', 'FPR-0.01', 'FPR-0.001', 'FPR-0.002', 'FPR-0.005']
fig, axs = plt.subplots(5, 1, figsize = (10,10))
mask = pred_nf_counts['IsX4']
for col, ax in zip(tcols, axs.flatten()):
    nbins = min(10, len(pred_nf_counts[col]))
    
    _, obins, _ = ax.hist([pred_nf_counts[col][mask],
                           pred_nf_counts[col][~mask]], 
                          bins = nbins,
                          label = ['X4', 'R5'])
    if ax.is_first_row():
        ax.legend()
    ax.set_title(col)
plt.tight_layout()

# <codecell>

env_seqs = pseqdf.ix['C'][['gp120', 'gp41']].dropna(how = 'all').fillna('-').apply(lambda x: ''.join(x), axis = 1)
ltr_seqs = pseqdf.ix['C']['ltr'].dropna()
vpr_seqs = pseqdf.ix['C']['vpr'].dropna()
tat_seqs = pseqdf.ix['C']['tat'].dropna()
nef_seqs = pseqdf.ix['C']['nef'].dropna()

subc_seqs = pd.DataFrame({
                          'LTR':ltr_seqs,
                          'Env':env_seqs,
                          'Vpr':vpr_seqs,
                          'Tat':tat_seqs,
                          'Nef':nef_seqs
                          })
subc_seqs, subc_pred_x4 = subc_seqs.align(x4_pred.ix['C'], axis = 0, join = 'inner')

# <codecell>

env_seqs = pseqdf.ix['B'][['gp120', 'gp41']].dropna(how = 'all').fillna('-').apply(lambda x: ''.join(x), axis = 1)
ltr_seqs = pseqdf.ix['B']['ltr'].dropna()
vpr_seqs = pseqdf.ix['B']['vpr'].dropna()
tat_seqs = pseqdf.ix['B']['tat'].dropna()
nef_seqs = pseqdf.ix['B']['nef'].dropna()

subb_seqs = pd.DataFrame({
                          'LTR':ltr_seqs,
                          'Env':env_seqs,
                          'Vpr':vpr_seqs,
                          'Tat':tat_seqs,
                          'Nef':nef_seqs
                          })
subb_seqs, subb_pred_x4 = subb_seqs.align(x4_pred.ix['B'], axis = 0, join = 'inner')

# <codecell>

conc_seqs = {}

conc_seqs['LTR'] = """tggaagggttaatttactctaagaaaaggcaagagatccttgatttgtgg
gtctatcacacaCaAGgctactTcCCtGAtTGGCAaaacTacACACCgGG
aCCaGGggtcAgatacCCacTgACctttGGaTGGtgcTtcAAgctaGTaC
CAgTtgacCCaaggGaagtaGAagAggccaacgaaggaGAaaacaaCtGt
tTgcTaCAcCCtatgagcCagCatGGaatgGAggAtgaacacagagAAgt
atTaaagTGGaagtttgacagtcaccTaGcacgcagacacatggcCcGcg
aGctacATCcGGAgTatTAcAaagACTGCTGacACagaaGGgACTtTccg
ctgggactttccactgGGGcgttccaggaggtgtggtctGGGCGGgActG
GGAGTGGtcaaCCCtCaGatGctgCATATAAGCagcTGCTTTtcgcctgt
actgggtctctctaggtagaccagatctgagcctgggagctctctggcta
tctagggaacccactgcttaagcctcaataaagcttgccttgagtgctct
aagtagtgtgtgccctctgttttgactctggtaactagagatccctcaga
cccttttggtagtgaggaaatctctagca""".replace('\n', '').upper()

conc_seqs['Vpr'] = """MEQAPEDQGPQREPYNEWTLELLEELKQEAVRHFPRPWLHSLGQYIYETY
GDTWTGVEAIIRILQQLLFIHFRIGCQHSRIGILRQRRARNGASRS""".replace('\n', '').upper()

conc_seqs['Tat'] = """MEPVDPNLEPWNHPGSQPKTACNKCYCKHCSYHCLVCFQTKGLGISYGRK
KRRQRRSAPPSSEDHQNPISKQPLPQTRGDPTGSEESKKKVESKTETDPFD""".replace('\n', '').upper()

conc_seqs['Nef'] = """MGGKWSKSSIVGWPAVRERIRRTEPAAEGVGAASQDLDKHGAL
TSSNTATNNADCAWLEAQEEEEEVGFPVRPQVPLRPMTYKAAFDLSFFLK
EKGGLEGLIYSKKRQEILDLWVYHTQGYFPDWQNYTPGPGVRYPLTFGWC
FKLVPVDPREVEEANEGENNCLLHPMSQHGMEDEDREVLKWKFDSHLARR
HMARELHPEYYKDC""".replace('\n', '').upper()

conc_seqs['Env'] = """MRVRGILRNCQQWWIWGILGFWMLMICNVVGNLWVTVYYGVPVWKEAKTT
LFCASDAKAYEKEVHNVWATHACVPTDPNPQEIVLENVTENFNMWKNDMV
DQMHEDIISLWDQSLKPCVKLTPLCVTLNCTNATNATNTM
GEIKNCSFNITTELRDKKQKVYALFYRLDIVPLNENNSY
RLINCNTSAITQACPKVSFDPIPIHYCAPAGYAILKCNNKTFNGTGPCNN
VSTVQCTHGIKPVVSTQLLLNGSLAEEEIIIRSENLTNNAKTIIVHLNES
VEIVCTRPNNNTRKSIRIGPGQTFYATGDIIGDIRQAHCNISEDKWNKTL
QKVSKKLKEHFPNKTIKFEPSSGGDLEITTHSFNCRGEFFYCNTSKLFNS
TYNSTNSTITLPCRIKQIINMWQEVGRAMYAPPIAGNIT
CKSNITGLLLTRDGGKNNTETFRPGGGDMRDNWRSELYKYKVVEIKP
LGIAPTKAKRRVVEREKRAVGIGAVFLGFLGAAGSTMGAASITLTVQARQ
LLSGIVQQQSNLLRAIEAQQHMLQLTVWGIKQLQTRVLAIERYLKDQQLL
GIWGCSGKLICTTAVPWNSSWSNKSQEDIWDNMTWMQWDREISNYTDTIY
RLLEDSQNQQEKNEKDLLALDSWKNLWNWFDITNWLWYIKIFIMIVGGLI
GLRIIFAVLSIVNRVRQGYSPLSFQTLTPNPRGPDRLGRIEEEGGEQDR
DRSIRLVSGFLALAWDDLRSLCLFSYHRLRDFILIAARAVELLGRSSLRG
LQRGWEALKYLGSLVQYWGLELKKSAISLLDTIAIAVAEGTDRIIELIQR
ICRAIRNIPRRIRQGFEAALQ""".replace('\n', '').upper()





# <codecell>

conb_seqs = {}

conb_seqs['LTR'] = """tGGAaGGGcTaaTttacTCcCaaaaaAgacAagAtATcCTTGAtcTGTGG
gtctAccAcAcaCAaGGCTacTTCCCtGAtTgGCAgAAcTACACAccAGG
gccaGGgatcAGaTatCCacTgaccTTTGGaTGgTGcTtcAAgcTAGTAC
CAgTtgAgCcAGagaaggTagAagagGccAatgaaggagagaacaacagc
tTGtTaCAcCCtatgagCCtgCATGGgatgGAgGAcccgGAgaaagAAGt
gtTagtgTGGAagTttGACAgccgccTaGcatttcatCAcatggCccgaG
AgctgcATCCggAgTactacaaggActGcTGACatcgagctttctacaaG
GGACTTTCCgCtgGGGACTTTccagggagGcGtggcctGGgcgggaCtgg
GGAgtggCgagCCCtcAGAtgcTgCATATAAGCAGCtGCttTttGccTGT
ACtGGgTCTCTCTggttaGaCCAGATCtGAGCctGGGAGcTCtctggcta
actagggaacccactgcttaagcctcaataaagcttgccttgagtgcttc
aagtagtgtgtgcccgtctgttgtgtgactctggtaactagagatccctc
agacccttttagtcagtgtggaaaatctctagca""".replace('\n', '').upper().replace('-', '')

conb_seqs['Vpr'] = """MEQAPEDQGPQREPYNEWTLELLEELKNEAVRHFPRIWLHGLGQHIYETY
GDTWAGVEAIIRILQQLLFIHFRIGCQHSRIGIT-RQRRARNGASRS""".replace('\n', '').upper().replace('-', '')

conb_seqs['Tat'] = """MEPVDPRLEPWKHPGSQPKTACTNCYCKKCCFHCQVCFITKGLGISYGRK
KRRQRRRAPQDSQTHQVSLSKQPASQ-PRGDPTGPKESKKKVERETETDP
VD""".replace('\n', '').upper().replace('-', '')

conb_seqs['Nef'] = """MGGKWSKRSVVGWPTVRERMR-----RAEPAA--DGVGAVSRDLEKHGAI
TSSNTAANNADCAWLEAQEEE-EVGFPVRPQVPLRPMTYKGALDLSHFLK
EKGGLEGLIYSQKRQDILDLWVYHTQGYFPDWQNYTPGPGIRYPLTFGWC
FKLVPVEPEKVEEANEGENNSLLHPMSLHGMDDPEREVLVWKFDSRLAFH
HMARELHPEYYKDC""".replace('\n', '').upper().replace('-', '')

conb_seqs['Env'] = """MRVMGTQRNYQHLWRWGILILGMLIMCKAT-DLWVTVYYGVPVWKDADTT
LFCASDAKAYDTEVHNVWATHACVPTDPNPQEVNLENVTEDFNMWKNNMV
EQMHEDIISLWDQSLKPCVKLTPLCVTLNCSNA--NTTN--NSTM-----
----EEIKNCSYNITTELRDKTQKVYSLFYKLDVVQLDES--NKSEYY-Y
RLINCNTSAITQACPKVSFEPIPIHYCAPAGFAILKCKDPRFNGTGSCNN
VSSVQCTHGIKPVASTQLLLNGSLAEGKVMIRSENITNNAKNIIVQFNKP
VPITCIRPNNNTRKSIRFGPGQAFYT-NDIIGDIRQAHCNINKTKWNATL
QKVAEQLREHFPNKTIIFTNSSGGDLEITTHSFNCGGEFFYCNTTGLFNS
TWK-NGTT---NNTEQM--ITLPCRIKQIINMWQRVGRAMYAPPIAGVIK
CTSNITGIILTRDGGNNET---ETFRPGGGDMRDNWRSELYKYKVVKIEP
LGVAPTRAKRRVVEREKRAVGMGAVFLGFLGAAGSTMGAASITLTVQARQ
LLSGIVQQQSNLLKAIEAQQHLLKLTVWGIKQLQARVLALERYLQDQQLL
GIWGCSGKLICATTVPWNSSWSNKTQEEIWNNMTWLQWDKEISNYTNIIY
KLLEESQNQQEKNEQDLLALDKWANLWNWFNITNWLWYIRIFIMIVGGLI
GLRIVIAIISVVNRVRQGYSPLSFQIPTP-NPEGLDRPGRIEEGGGEQGR
DRSIRLVSGFLALAWDDLRSLCLFSYHRLRDCILIAARTVELLGHSSLKG
LRLGWEGLKYLWNLLLYWGRELKNSAISLLDTIAVAVAEWTDRVIEIGQR
ACRAILNIPRRIRQGFERALL""".replace('\n', '').upper().replace('-', '')

# <codecell>

aligns = []
for col in subc_seqs.columns:
    
    seqs = subc_seqs[col].dropna().to_dict().items()
    for num, (gi, align) in enumerate(GeneralSeqTools.seq_align_to_ref(seqs, conc_seqs[col])):
        if (num == 0) or (num == 100) or (num == 10) or (num % 1000 == 0):
            print col, num
            
        aligns.append({
                       'Prot':col,
                       'ID':gi,
                       'Align':align
                       })

# <codecell>

baligns = []
for col in subb_seqs.columns:
    
    seqs = subb_seqs[col].dropna().to_dict().items()
    for num, (gi, align) in enumerate(GeneralSeqTools.seq_align_to_ref(seqs, conb_seqs[col])):
        if (num == 0) or (num == 100) or (num == 10) or (num % 1000 == 0):
            print col, num
            
        baligns.append({
                       'Prot':col,
                       'ID':gi,
                       'Align':align
                       })

# <codecell>

subc_align = pd.pivot_table(pd.DataFrame(aligns),
                               rows = 'ID',
                               cols = 'Prot',
                               values = 'Align',
                               aggfunc = 'first')

# <codecell>

subb_align = pd.pivot_table(pd.DataFrame(baligns),
                               rows = 'ID',
                               cols = 'Prot',
                               values = 'Align',
                               aggfunc = 'first')

# <codecell>

tmp_seqs = []
for col in subc_align.columns:
    res = subc_align[col].fillna('-').apply(lambda x: pd.Series(list(x)))
    res.columns = pd.MultiIndex.from_tuples([(col, pos+1) for pos in res.columns],
                                          names = ['Prot', 'OCol'])
    tmp_seqs.append(res.copy())
    
subc_cols = pd.concat(tmp_seqs, axis = 1).replace('-', np.nan).replace('X', np.nan)
subc_cols, _ = subc_cols.align(subc_seqs, axis = 0, join = 'right')

# <codecell>

tmp_seqs = []
for col in subb_align.columns:
    res = subb_align[col].fillna('-').apply(lambda x: pd.Series(list(x)))
    res.columns = pd.MultiIndex.from_tuples([(col, pos+1) for pos in res.columns],
                                          names = ['Prot', 'OCol'])
    tmp_seqs.append(res.copy())
    
subb_cols = pd.concat(tmp_seqs, axis = 1).replace('-', np.nan).replace('X', np.nan)
subb_cols, _ = subb_cols.align(subb_seqs, axis = 0, join = 'right')

# <codecell>

from Bio.Alphabet import generic_dna, generic_protein
from itertools import product, islice, imap
from collections import deque
import os, os.path
import TreeingTools

alpha = generic_dna
def append_seq(ser):
    return ''.join(ser.fillna('-'))


def rolling_tree_apply(tup):
    
    group_series, seq_series, kwargs = tup
    
    fname = '/home/will/SubCData/Trees/Tree-%(sub)s-%(Prot)s-%(Start)i-%(WinSize)i.newick' % kwargs
    if os.path.exists(fname):
        return True
    
    
    alpha = generic_dna if kwargs['Prot'] == 'LTR' else generic_protein
    
    seq_series = seq_series.dropna(thresh = 5)
    
    vseq, vgroup = seq_series.align(group_series.dropna(), join = 'inner', axis = 0)
    
    nseq_ser = vseq.apply(append_seq, axis = 1)
    nseqs = sorted(nseq_ser.to_dict().items())
    
    trop_dict = vgroup.to_dict()
    #print nseqs
    #try:
    #    tree, dmat = TreeingTools.phylip_tree_collapse_unique(nseqs, alphabet=alpha, use_fast=True)
    #except:
    #    return False
    #print 'treeing', fname
    tree = TreeingTools.run_FastTree(nseqs, alphabet=alpha, uniq_seqs=True)
    
    with open(fname, 'w') as handle:
        tree.write(handle, schema='newick')
    return True
    
    
    try:
        tree, dmat = TreeingTools.phylip_tree_collapse_unique(nseqs, alphabet=alpha, use_fast=True)
        benj_res = TreeingTools.check_distance_pvals(dmat, trop_dict, nreps = 50)
    except:
        return kwargs
    
    benj_res.update(kwargs)
    try:
        out = TreeingTools.evaluate_association_index(tree, trop_dict)
        benj_res['AI'], benj_res['AI-pval'], benj_res['AI-null'] = out
    except:
        benj_res['AI'], benj_res['AI-pval'], benj_res['AI-null'] = (None, None, None)
        
    return benj_res


def yield_regions(seq_data, trop_data, windows, prots, sub):
    
    for win, prot in product(windows, prots):
        print win, prot
        tmp_df = seq_data[prot]
        for start in range(1, len(tmp_df.columns)-win):
            seqs = tmp_df[range(start, start+win)]
            tdict = {
                     'WinSize':win,
                     'Start':start,
                     'Prot':prot,
                     'sub':sub
                     }
            yield trop_data, seqs, tdict

# <codecell>

raise KeyboardInterrupt

# <codecell>

iterable = yield_regions(subc_cols, subc_pred_x4, [35, 10, 15], ['Env', 'Vpr', 'Tat', 'LTR', 'Nef'], 'C')
with ProcessPoolExecutor(max_workers = 30) as pool:
    for num, res in enumerate(pool.map(rolling_tree_apply, iterable)):
        if (num == 0) or (num == 100) or (num == 1000) or (num % 10000 == 0):
            print num, res

# <codecell>

iterable = yield_regions(subb_cols, subb_pred_x4, [35, 10, 15], ['Env', 'Vpr', 'Tat', 'LTR', 'Nef'], 'B')
with ProcessPoolExecutor(max_workers = 30) as pool:
    for num, res in enumerate(pool.map(rolling_tree_apply, iterable)):
        if (num == 0) or (num == 100) or (num == 1000) or (num % 10000 == 0):
            print num, res

# <codecell>

from collections import defaultdict
from itertools import product, combinations
from scipy.stats import ttest_ind
from random import shuffle

def check_root(tree, trop_dict, nreps):
    groups = defaultdict(list)
    missing = 0
    for node in tree.leaf_nodes():
        try:
            groups[trop_dict[node.taxon.label]].append(node.distance_from_root())
        except KeyError:
            missing += 1
    _, pval = ttest_ind(groups[True], groups[False])
    
    min_size = min(len(groups[True]), len(groups[False]))-1
    pvals = []
    for n in range(nreps):
        shuffle(groups[True])
        shuffle(groups[False])
        _, tpval = ttest_ind(groups[True][:min_size], groups[False][:min_size])
        pvals.append(tpval)
    
    return pval, np.mean(pvals)

# <codecell>

import dendropy
tree_files = sorted(glob.glob('/home/will/SubCData/Trees/Tree-C*.newick'))
#fname = '/home/will/SubCData/Trees/Tree-%(Prot)s-%(Start)i-%(WinSize)i.newick' % kwargs
results = []
trop_dict = subc_pred_x4.to_dict()
for num, f in enumerate(tree_files):
    if (num == 0) or (num == 100) or (num % 500) == 0:
        print num, f
    _, sub, prot, start, winsize = f.rsplit(os.sep, 1)[-1].rsplit('.',1)[0].split('-')
    tree = dendropy.Tree.get_from_path(f, 'newick')
    try:
        pval, adj_pval = check_root(tree, trop_dict, 20)
    except ZeroDivisionError:
        pval, adj_pval = (np.nan, np.nan)
    
    results.append({
                    'Prot':prot,
                    'Start': float(start),
                    'Winsize': float(winsize),
                    'Pvalue': pval,
                    'AdjPval': adj_pval})
                    

# <codecell>

res_df = pd.DataFrame(results).groupby(['Prot', 'Winsize', 'Start']).first()
res_df['LogAdj'] = -np.log10(res_df['AdjPval'])
res_df

# <codecell>

t = res_df.copy().reset_index()
t['Midpoint'] = t['Start']+t['Winsize']/2

t.to_excel('/home/will/Dropbox/Wigdahl HIV Lab/SubCAnalysis/position_pvals.xlsx')

# <codecell>

prots = results['Prot'].unique()
wins = results['Winsize'].unique()
print prots, wins

# <codecell>

fig, axs = plt.subplots(5,1, figsize = (10, 15))
prots = ['LTR', 'Vpr', 'Nef', 'Tat', 'Env']
wins = [10, 15, 35]
colors = 'rgb'
for ax, prot in zip(axs.flatten(), prots):
    
    for win, c in zip(wins, colors):
        tmp = res_df.ix[prot].ix[win]
        ax.plot(tmp.index, tmp['LogAdj'], c, alpha = 0.7, label = 'Winsize-%i' % win)
    ax.set_title(prot)
    ax.hlines(2, 0, len(conc_seqs[prot]))
    ax.set_xlim(0, len(conc_seqs[prot]))
    if prot == 'Env':
        ax.set_yscale('log')
    if ax.is_first_row():
        ax.legend()
fig.tight_layout()
plt.savefig('/home/will/Dropbox/Wigdahl HIV Lab/SubCAnalysis/tree_analysis_raw.png')

# <codecell>

fig, axs = plt.subplots(5,1, figsize = (10, 15))
prots = ['LTR', 'Vpr', 'Nef', 'Tat', 'Env']
wins = [10, 15, 35]
colors = 'rgb'
for ax, prot in zip(axs.flatten(), prots):
    
    for win, c in zip(wins, colors):
        tmp = res_df.ix[prot].ix[win]
        ax.plot(tmp.index, pd.rolling_max(tmp['LogAdj'], win, min_periods = 1), 
                c, alpha = 0.7, label = 'Winsize-%i' % win)
    ax.set_title(prot)
    ax.hlines(2, 0, len(conc_seqs[prot]))
    ax.set_xlim(0, len(conc_seqs[prot]))
    if prot == 'Env':
        ax.set_yscale('log')
    if ax.is_first_row():
        ax.legend()
fig.tight_layout()
plt.savefig('/home/will/Dropbox/Wigdahl HIV Lab/SubCAnalysis/tree_analysis_rectified.png')

# <codecell>

typs = ['raw', 'rectified']
base_fname = '/home/will/Dropbox/Wigdahl HIV Lab/SubCAnalysis/TreeDiffFigs/'
#
for prot, win, typ in product(['gp120', 'gp41']+prots, wins, typs):
    fig, ax = plt.subplots(1,1, figsize = (10, 5))
    if prot == 'gp41':
        tmp = res_df.ix['Env'].ix[win].ix[511:].reset_index(drop=True).copy()
    elif prot == 'gp120':
        tmp = res_df.ix['Env'].ix[win].ix[:511].copy()
    else:
        tmp = res_df.ix[prot].ix[win].copy()
    
    reject, adjpval, _, _ = multipletests(tmp['Pvalue'].dropna(), method = 'bonferroni', alpha = 0.01)
    pval_threshold = adjpval[reject].max()
    
    if typ == 'rectified':
        tmp['LogAdj'] = pd.rolling_max(tmp['LogAdj'], win, min_periods = 1)
    ax.plot(np.array(tmp.index)+(win/2), tmp['LogAdj'])
    
    
    if prot == 'gp41':
        ax.set_xlim(0, len(conc_seqs['Env'])-511)
        
    elif prot == 'gp120':
        ax.set_xlim(0, 511)
    else:
        ax.set_xlim(0, len(conc_seqs[prot]))
    ax.set_ylim(0, ax.get_ylim()[1])
    ax.set_title('%s win-%i %s' % (prot, win, typ))
    ax.hlines(-np.log10(pval_threshold), *ax.get_xlim())
    
    tdict = {
             'prot':prot,
             'win':win,
             'typ':typ
             }
    fname = base_fname+'%(prot)s-%(win)i-%(typ)s.png' % tdict
    
    plt.savefig(fname, dpi = 300)
    plt.close()
    
    

# <codecell>

for prot, typ in product(['gp120', 'gp41']+prots, typs):
    fig, ax = plt.subplots(1,1, figsize = (10, 5))
    
    for win in wins:
        if prot == 'gp41':
            tmp = res_df.ix['Env'].ix[win].ix[511:].reset_index(drop=True).copy()
        elif prot == 'gp120':
            tmp = res_df.ix['Env'].ix[win].ix[:511].copy()
        else:
            tmp = res_df.ix[prot].ix[win].copy()
    
        reject, adjpval, _, _ = multipletests(tmp['Pvalue'].dropna(), method = 'bonferroni', alpha = 0.01)
        pval_threshold = adjpval[reject].max()
    
        if typ == 'rectified':
            tmp['LogAdj'] = pd.rolling_max(tmp['LogAdj'], win, min_periods = 5)
        ax.plot(np.array(tmp.index)+(win/2), tmp['LogAdj'])
    
    
    if prot == 'gp41':
        ax.set_xlim(0, len(conc_seqs['Env'])-511)
        
    elif prot == 'gp120':
        ax.set_xlim(0, 511)
    else:
        ax.set_xlim(0, len(conc_seqs[prot]))
    ax.set_ylim(0, ax.get_ylim()[1])
    ax.set_title('%s -win-%i' % (prot, win))
    ax.hlines(-np.log10(pval_threshold), *ax.get_xlim())
    
    tdict = {
             'prot':prot,
             'win':'ALL',
             'typ':typ
             }
    fname = base_fname+'%(prot)s-%(win)s-%(typ)s.png' % tdict
    
    plt.savefig(fname, dpi = 300)
    plt.close()

# <codecell>

raise AssertionError

# <codecell>

ld = []
for f in glob.glob('/home/will/SubCData/LANLRes/LANLResults_*.txt'):
    ld.append(pd.read_csv(f, sep='\t', index_col=0))
    
lanl_data = pd.concat(ld, axis = 0, ignore_index=True)

# <codecell>

from types import StringType
def norm_id(inser):
    
    cols = ['Accession', 'GI number']
    for col in cols:
        if type(inser[col]) == StringType:
            return inser[col].split('.')[0]
    print inser
    raise KeyboardInterrupt
    

lanl_data['GBid'] = lanl_data.apply(norm_id, axis = 1)
agg_lanl_data = lanl_data.groupby('GBid').first()

# <codecell>

cols, ex_data = subc_cols.align(agg_lanl_data, join = 'left', axis = 0)

# <codecell>

agg_lanl_data

# <codecell>

seqdf['HasSeq'] = 1.0
seq_counts = pd.pivot_table(seqdf,
                            rows = 'ID',
                            cols = 'Prot',
                            values = 'HasSeq',
                            aggfunc = 'sum')
ex_data, seq_count_aln = agg_lanl_data.align(seq_counts, axis = 0, join = 'inner')

# <codecell>

has_data = ex_data.copy()
for col in ex_data.columns:
    has_data[col] = ex_data[col].notnull()

#pd.crosstab(seq_counts

# <codecell>

cols = ['CD4 count', 'CD8 count', 'Coreceptor', 
        'Days from Infection', 'Days from Seroconversion',
        'Patient Age', 'Viral load', 'Patient Sex', 'Patient Sex']
counts = []
for prot1, prot2, col in product(seq_count_aln.columns, seq_count_aln.columns, cols):
    
    m1 = seq_count_aln[prot1]
    m2 = seq_count_aln[prot2]
    hd = has_data[col]
    
    mask = m1 * m2 * hd
    
    counts.append({
                   'Prot1':prot1,
                   'Prot2':prot2,
                   'Measure': col,
                   'Count':mask.sum()
                   
                   })

# <codecell>

type_counts = pd.pivot_table(pd.DataFrame(counts),
                             rows = ['Measure', 'Prot1'],
                             cols = 'Prot2',
                             values = 'Count',
                             aggfunc = 'sum',
                             margins = True)
type_counts.to_excel('/home/will/Dropbox/Wigdahl HIV Lab/SubCAnalysis/data_counts.xlsx')

# <codecell>

type_counts.head()

# <codecell>

base_fname

# <codecell>


