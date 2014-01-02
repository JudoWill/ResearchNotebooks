# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import os, os.path
from matplotlib import pyplot as plt
from pylab import get_cmap
import SimpleCV as cv
from glob import glob

# <codecell>

def show_img(img, ax = None):
    if ax is not None:
        plt.sca(ax)
    nimg = img.getNumpy()
    return plt.imshow(nimg, aspect='equal')

# <codecell>

path = '/home/will/Dropbox/burnimages/*.jpg'
norm_files = sorted(f for f in glob(path) if '-e' not in f)
masked_files = sorted(f for f in glob(path) if '-e' in f)

fig, axs = plt.subplots(6,6, figsize = (10,10))
for f, ax in zip(norm_files, axs.flatten()):
    img = cv.Image(f)
    show_img(img, ax = ax)
    ax.set_xticks([])
    ax.set_yticks([])
fig.tight_layout()

# <codecell>

from itertools import islice, izip_longest
from dateutil.parser import parse
def make_wound_mask(norm_img, green_img, color,
                    minsize = None,
                    maxsize = None):
    
    wmask = green_img.hueDistance(color).invert().threshold(200)
    blobs = norm_img.findBlobsFromMask(wmask, 
                                       minsize = minsize,
                                       maxsize = maxsize)
    
    return wmask, blobs

fig, axs = plt.subplots(6,6, figsize = (10,10))
results = []
for fname, mf, of, ax in izip_longest(norm_files, masked_files, norm_files, axs.flatten()):
    mask_img = cv.Image(mf)
    norm_img = cv.Image(of)
    
    dt = parse(fname.rsplit(os.sep,1)[1].replace('.jpg', '').replace('.',':'))

    wound_mask, wound_blobs = make_wound_mask(norm_img, mask_img, cv.Color.GREEN,
                                              minsize = 1000)
    dime_mask, dime_blobs = make_wound_mask(norm_img, mask_img, cv.Color.BLUE,
                                            minsize = 500)
    layer = cv.DrawingLayer((norm_img.width, norm_img.height))
    wound_blobs[-1].drawHull(color=cv.Color.BLUE, width = 100, layer = layer)
    dime_blobs[-1].drawHull(color=cv.Color.RED, width = 100, layer = layer)
    
    norm_img.addDrawingLayer(layer)
    fnorm = norm_img.applyLayers()
    
    ratio = wound_blobs[-1].area()/dime_blobs[-1].area()
    results.append((dt, ratio))
    if ax is not None:
        show_img(fnorm, ax = ax)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(ratio)
        
fig.tight_layout()

# <codecell>

import pandas as pd

res_df = pd.DataFrame(sorted(results), columns = ['SampleTime', 'Ratio'])
dime_diameter = 18 #mm
dime_area = 3.141*(dime_diameter/2)**2
res_df['Area-mm2'] = dime_area*res_df['Ratio']
res_df.set_index('SampleTime', inplace=True)
res_df

# <codecell>

res_df['Area-mm2'].plot()
out = pd.ewma(res_df['Area-mm2'], freq='d', span = 1)
out.plot(lw = 10, alpha = 0.7)
plt.ylabel('Wound-Area-mm^2')

# <codecell>


