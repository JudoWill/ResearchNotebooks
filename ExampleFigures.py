# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

from __future__ import division
import os, os.path
from pandas import *
import numpy as np
import csv
import matplotlib.pyplot as plt
from itertools import product, imap, islice
from patsy import dmatrices
from patsy.contrasts import Treatment
import statsmodels.api as sm
import sys

os.chdir('/home/will/HIVSystemsBio/NewCytokineAnalysis/')

# <codecell>

points = np.abs((5*np.random.randn(3, 50)+np.tile(np.arange(1,51), (3, 1))).transpose())

from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.axis3d as axis3d
ncolor = (103/255, 108/255, 124/255, 0.75)
# New axis settings
custom_AXINFO = {
    'x': {'i': 0, 'tickdir': 1, 'juggled': (1, 0, 2),
          'color': ncolor},
    'y': {'i': 1, 'tickdir': 0, 'juggled': (0, 1, 2),
          'color': ncolor},
    'z': {'i': 2, 'tickdir': 0, 'juggled': (0, 2, 1),
          'color': ncolor},}

class custom_XAxis(axis3d.Axis):
    _AXINFO = custom_AXINFO

class custom_YAxis(axis3d.Axis):
    _AXINFO = custom_AXINFO

class custom_ZAxis(axis3d.Axis):
    _AXINFO = custom_AXINFO

class custom_Axes3D(Axes3D):
    def _init_axis(self):
        '''Init 3D axes; overrides creation of regular X/Y axes'''
        self.w_xaxis = custom_XAxis('x', self.xy_viewLim.intervalx,
                                    self.xy_dataLim.intervalx, self)
        self.xaxis = self.w_xaxis
        self.w_yaxis = custom_YAxis('y', self.xy_viewLim.intervaly,
                            self.xy_dataLim.intervaly, self)
        self.yaxis = self.w_yaxis
        self.w_zaxis = custom_ZAxis('z', self.zz_viewLim.intervalx,
                            self.zz_dataLim.intervalx, self)
        self.zaxis = self.w_zaxis

        for ax in self.xaxis, self.yaxis, self.zaxis:
            ax.init3d()


# <codecell>

from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure(figsize = (10,10))
ax = custom_Axes3D(fig)
ax.scatter(points[:,0], points[:,1], points[:,2], s = 45, c='red')
ax.view_init(elev=0., azim=0)
ax.set_xticks([]);
ax.set_ylim([0, 60])
ax.set_zlim([0, 60])
ax.set_xlim([0, 60])
ax.set_zlabel('Cytokine')
ax.set_ylabel('Parameter')
ax.grid(linewidth=20)
plt.savefig("figures/example_images/img-1.png", dpi = 300)

# <codecell>

from scipy.stats import linregress
m, b, _, _, _ = linregress(points[:,1], points[:,2])
reg_x = np.zeros((51,1))
reg_y = np.arange(51)
reg_z = m*reg_y+b

fig = plt.figure(figsize = (10,10))
ax = custom_Axes3D(fig)
ax.scatter(points[:,0], points[:,1], points[:,2], s = 45, c='red')
ax.plot3D(reg_x, reg_y, reg_z)
ax.view_init(elev=0., azim=0)
ax.set_xticks([]);
ax.set_ylim([0, 60])
ax.set_zlim([0, 60])
ax.set_xlim([0, 60])
ax.set_zlabel('Cytokine')
ax.set_ylabel('Parameter')
plt.savefig("figures/example_images/img-2.png", dpi = 300)

# <codecell>

fig = plt.figure(figsize = (10,10))
ax = custom_Axes3D(fig)
ax.scatter(points[:,0], points[:,1], points[:,2], s = 45, c='red')
#ax.contourf(points[:,0], points[:,1], points[:,2])
ax.plot3D(reg_x, reg_y, reg_z)
#ax.view_init(elev=20, azim=20)
ax.set_ylim([0, 60])
ax.set_zlim([0, 60])
ax.set_xlim([0, 60])
ax.set_zlabel('Cytokine')
ax.set_ylabel('Parameter-1')
ax.set_xlabel('Parameter-2')
plt.savefig("figures/example_images/img-3.png", dpi = 300)

# <codecell>

nreg_x = np.arange(51)
nreg_y = np.arange(51)
nreg_z = m*(nreg_x+nreg_y)/2+b


fig = plt.figure(figsize = (10,10))
ax = custom_Axes3D(fig)
ax.scatter(points[:,0], points[:,1], points[:,2], s = 45, c='red')
ax.plot3D(nreg_x, nreg_y, nreg_z)
#ax.view_init(elev=20, azim=20)
ax.set_ylim([0, 60])
ax.set_zlim([0, 60])
ax.set_xlim([0, 60])
ax.set_zlabel('Cytokine')
ax.set_ylabel('Parameter-1')
ax.set_xlabel('Parameter-2')
plt.savefig("figures/example_images/img-4.png", dpi = 300)

# <codecell>

nvals = 10*np.abs(np.random.randn(10,5))
df = DataFrame(nvals)
df[3] *= 10
df.ix[2] *= 6

# <codecell>

plt.figure(figsize = (10,10))
ax = plt.subplot(111)
ncols = 5
nrows = 10
for col in range(0,ncols):
    for row in range(0, nrows):
        val = '%.02f' % (df.values[row, col])
        ax.text(col/ncols+(1/ncols)/2, row/nrows+0.03, val,
                 transform = ax.transAxes,
                 fontsize = 30,
                 horizontalalignment='center')
        ax.add_patch(Rectangle((col/ncols, row/nrows), 
                                height=1/nrows, 
                                width=1/ncols,
                                facecolor = 'white'))
ax.set_xticks([])
ax.set_yticks([])
plt.title('Raw data')
plt.savefig('figures/example_images/qimg-1.png', dpi = 300)

# <codecell>

plt.figure(figsize = (10,10))
ax = plt.subplot(111)

ranks = df.rank()

for col in range(0,ncols):
    for row in range(0, nrows):
        val = '%.02f' % (df.values[row, col])
        ax.text(col/ncols+(1/ncols)/2, row/nrows+0.03, val,
                 transform = ax.transAxes,
                 fontsize = 30,
                 horizontalalignment='center')
        ax.add_patch(Rectangle((col/ncols, row/nrows), 
                                height=1/nrows, 
                                width=1/ncols,
                                alpha=0.4,
                                facecolor = [0.5, ranks.values[row, col]/11, 1]))
plt.title('Ranked data')
plt.savefig('figures/example_images/qimg-2.png', dpi = 300)

# <codecell>

plt.figure(figsize = (10,10))
ax = plt.subplot(111)

ranks = df.rank()

for col in range(0,ncols):
    for row in range(0, nrows):
        val = '%.02f' % (df.values[row, col])
        ax.text(col/ncols+(1/ncols)/2, row/nrows+0.03, val,
                 transform = ax.transAxes,
                 fontsize = 30,
                 horizontalalignment='center')
        edgecolor = 'red' if ranks.values[row, col] == 1 else 'black'
        lw = 10 if ranks.values[row, col] == 1 else 2
        ax.add_patch(Rectangle((col/ncols, row/nrows), 
                                height=1/nrows, 
                                width=1/ncols,
                                alpha=0.4,
                                edgecolor=edgecolor,
                                lw = lw,
                                facecolor=[0.5, ranks.values[row, col]/11, 1]))
plt.title('Selecting smallest item in each column')
plt.savefig('figures/example_images/qimg-3.png', dpi = 300)

# <codecell>

plt.figure(figsize = (10,10))
ax = plt.subplot(111)

ndf = df.copy()
mask = ranks.values == 1
mval = np.mean(df.values[mask])
ndf.values[mask] = mval
for col in range(0,ncols):
    for row in range(0, nrows):
        val = '%.02f' % (ndf.values[row, col])
        ax.text(col/ncols+(1/ncols)/2, row/nrows+0.03, val,
                 transform = ax.transAxes,
                 fontsize = 30,
                 horizontalalignment='center')
        edgecolor = 'red' if ranks.values[row, col] == 1 else 'black'
        lw = 10 if ranks.values[row, col] == 1 else 2
        ax.add_patch(Rectangle((col/ncols, row/nrows), 
                                height=1/nrows, 
                                width=1/ncols,
                                alpha=0.4,
                                edgecolor=edgecolor,
                                lw = lw,
                                facecolor=[0.5, ranks.values[row, col]/11, 1]))
plt.title('Replacing with average value')
plt.savefig('figures/example_images/qimg-4.png', dpi = 300)                

# <codecell>

plt.figure(figsize = (10,10))
ax = plt.subplot(111)

ranks = df.rank()

for col in range(0,ncols):
    for row in range(0, nrows):
        val = '%.02f' % (ndf.values[row, col])
        ax.text(col/ncols+(1/ncols)/2, row/nrows+0.03, val,
                 transform = ax.transAxes,
                 fontsize = 30,
                 horizontalalignment='center')
        if ranks.values[row, col] >= 2:
            facecolor = [0.5, ranks.values[row, col]/11, 1]
        else:
            facecolor = 'grey'
        edgecolor = 'red' if ranks.values[row, col] == 2 else 'black'
        lw = 10 if ranks.values[row, col] == 2 else 2
        ax.add_patch(Rectangle((col/ncols, row/nrows), 
                                height=1/nrows, 
                                width=1/ncols,
                                alpha=0.4,
                                edgecolor=edgecolor,
                                lw = lw,
                                facecolor=facecolor))
plt.title('Pick second lowest value')
plt.savefig('figures/example_images/qimg-5.png', dpi = 300)

# <codecell>

plt.figure(figsize = (10,10))
ax = plt.subplot(111)

mask = ranks.values == 2
mval = np.mean(df.values[mask])
ndf.values[mask] = mval
for col in range(0,ncols):
    for row in range(0, nrows):
        val = '%.02f' % (ndf.values[row, col])
        ax.text(col/ncols+(1/ncols)/2, row/nrows+0.03, val,
                 transform = ax.transAxes,
                 fontsize = 30,
                 horizontalalignment='center')
        if ranks.values[row, col] >= 2:
            facecolor = [0.5, ranks.values[row, col]/11, 1]
        else:
            facecolor = 'grey'
        edgecolor = 'red' if ranks.values[row, col] == 2 else 'black'
        lw = 10 if ranks.values[row, col] == 2 else 2
        ax.add_patch(Rectangle((col/ncols, row/nrows), 
                                height=1/nrows, 
                                width=1/ncols,
                                alpha=0.4,
                                edgecolor=edgecolor,
                                lw = lw,
                                facecolor=facecolor))
plt.title('Replaced with mean value')
plt.savefig('figures/example_images/qimg-6.png', dpi = 300)                

# <codecell>

for num in range(1, 11):
    mask = ranks.values == num
    mval = np.mean(df.values[mask])
    ndf.values[mask] = mval

plt.figure(figsize = (10,10))
ax = plt.subplot(111)
for col in range(0,ncols):
    for row in range(0, nrows):
        val = '%.02f' % (ndf.values[row, col])
        ax.text(col/ncols+(1/ncols)/2, row/nrows+0.03, val,
                 transform = ax.transAxes,
                 fontsize = 30,
                 horizontalalignment='center')
        facecolor = 'grey'
        edgecolor = 'black'
        lw = 2
        ax.add_patch(Rectangle((col/ncols, row/nrows), 
                                height=1/nrows, 
                                width=1/ncols,
                                alpha=0.4,
                                edgecolor=edgecolor,
                                lw = lw,
                                facecolor=facecolor))
        
plt.title('completely normalized')
plt.savefig('figures/example_images/qimg-7.png', dpi = 300)                

# <codecell>


