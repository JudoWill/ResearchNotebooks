# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from __future__ import division
import SimpleCV as cv
import os
from random import randint
import numpy as np
import pandas as pd
disp = cv.Display(displaytype='notebook')
os.chdir('/home/will/SimpleCVExample/')

# <codecell>

print 2+2

# <codecell>

img.findCircle?

# <codecell>

from sklearn.metrics import mean_absolute_error
def gen_random_images(ncircles):
    base_img = cv.Image('empty_background.png').resize(100, 100)
    for _ in range(ncircles):
        pos = (randint(10, 40), randint(10, 40))
        base_img.dl().circle(pos, 10, cv.Color.RED, filled=True)
    nimg = base_img.applyLayers()
    mask = nimg.colorDistance(cv.Color.RED).invert().threshold(100)
    true_area = (mask.getNumpy()[:,:,0].flatten()>0).sum()
    return nimg, true_area


def predict_image(img, display=False, threshval = 50,
                  minsize = 10, 
                  maxsize = 0, threshblocksize = 5, 
                  threshconstant = 5, appx_level = 3):
    
    mask = img.colorDistance(cv.Color.RED).invert()
    #blobs = mask.findCircle()
    blobs = mask.findBlobs(threshval=threshval,
                           minsize=minsize,
                           maxsize=maxsize,
                           threshblocksize=threshblocksize,
                           threshconstant=threshconstant,
                           appx_level=appx_level)
    if display:
        if blobs is not None:
            blobs.draw()
        mask.save(disp)
            
    if blobs is None:
        return 0
    return blobs.area().sum()


def predict_batch(nreps, ncircles, **kwargs):
    
    true_areas = []
    guessed_areas = []
    for n in range(nreps):
        img, area = gen_random_images(ncircles)
        guessed_area = predict_image(img, **kwargs)
        true_areas.append(area)
        guessed_areas.append(guessed_area)
        
    return mean_absolute_error(true_areas, guessed_areas)

# <codecell>

img.findBlobs?

# <codecell>

img, area = gen_random_images(3)
img.save(disp)
print area

# <codecell>

predict_image(img, display=True,
              threshval=50)

# <codecell>

list(product(minsizes, threshvals))

# <codecell>

from itertools import product

minsizes = range(5, 100, 10)
#threshblocksizes = range(5, 20, 2)
#threshconstants = range(5, 20, 2)
#appx_levels = range(1, 11, 2)
threshvals = range(20, 200, 20)
circles = range(2, 5)

reps = 20


results = []

for circle in circles:
    print circle
    iterable = product(minsizes, threshvals)
    for minsize, threshval in iterable:
    
        mean_error = predict_batch(reps, circle, 
                                   threshval=threshval,
                                   minsize=minsize)
        results.append({
                    'circle': circle,
                    'minsize':minsize,
                    'threshval':threshval,
                    'mean_error': mean_error,
                    })


# <codecell>

results_df = pd.DataFrame(results)
pos = results_df['mean_error'].idxmin()
results_df.head(n=5)

# <codecell>

plt.figure(figsize = (15,15))
plt.scatter(results_df['mean_error'], results_df['threshval'])

# <codecell>

results_df.ix[pos]

# <codecell>


# <codecell>

results_df.ix[pos]

# <codecell>

_ = pd.scatter_matrix(results_df, figsize=(15,15))

# <codecell>

results_df['circle'].unique()

# <codecell>

print 100**2

# <codecell>


