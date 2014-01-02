# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import numpy as np
import urllib2 as urllib
from PIL import Image
from cStringIO import StringIO
from scipy.ndimage import label
from scipy.ndimage.measurements import find_objects
from random import shuffle
import csv
import re


# <codecell>

def prepend_if_needed(item, prepend):
    if not item.startswith(prepend):
        return prepend + item
    return item

def get_genes_by_pathway(pathway):
    url = 'http://rest.kegg.jp/link/genes/'+pathway
    data = urllib.urlopen(url)
    for row in csv.reader(data, delimiter = '\t'):
        yield row[1]

def color_pathway_by_genes(pathway, wanted_genes, unwanted_genes):
    url = 'http://www.kegg.jp/kegg-bin/show_pathway?'
    url += pathway.replace('path:', '') + '/default%3dwhite'
    for gene in wanted_genes:
        url += gene + '%09red,black/'
    for gene in unwanted_genes:
        url += gene + '%09white,black/'
        
    robj = re.compile('\<img src=\"(.+\.png)\" name=\"pathwayimage\" usemap=\"\#mapdata\" border=\"0\" /\>')
    web = urllib.urlopen(url).read()
    tmp = robj.findall(web)[0]
                      
    return 'http://www.kegg.jp/'+tmp
    
    

def get_img_url_from_genes(pathway, genes):
    
    pathway = prepend_if_needed(pathway, 'path:')
    pathway_genes = list(get_genes_by_pathway(pathway))
    gene_prepend = pathway_genes[0].split(':',1)[0]
    gene_set = set(prepend_if_needed(gene, gene_prepend + ':') for gene in genes)
    
    unwanted_genes = set(pathway_genes) - gene_set
    wanted_genes = gene_set & set(pathway_genes)
    url = color_pathway_by_genes(pathway, wanted_genes, unwanted_genes)
    return color_pathway_by_genes(pathway, wanted_genes, unwanted_genes)

def get_base_image_url(pathway):
    return get_img_url_from_genes(pathway, [])
    

# <codecell>



def get_numpy_from_url(url):
    img_file = urllib.urlopen(url)
    im = StringIO(img_file.read())
    resized_image = Image.open(im)
    return numpy.array(resized_image)

def get_colored_regions(imgdata):
    s = np.ones((3,3))
    nmask = np.dstack([imgdata[:,:,0]==255, imgdata[:,:,1]==0, imgdata[:,:,2]==0])
    labeled, num_items = label(np.all(nmask, axis = 2), s)
    return labeled, num_items

def make_zebra_stripe_inds(start, stop, num_stripes):
    tinds = np.linspace(start, stop, num_stripes+1)
    for i1, i2 in zip(tinds, tinds[1:]):
        yield i1, i2
        
def combine_images(img_tup):
    
    image_stack = np.dstack((img>0)*ind for ind, img in enumerate(img_tup,1))
    image_num_items = np.sum(image_stack>0, axis = 2)
    baseimage = np.zeros_like(image_num_items)
    s = np.ones((3,3))
    labeled, num_items = label(image_num_items>0, s)
    objects = find_objects(labeled)
    
    for sli in objects:
        img_present = [(img[sli]>1)*ind for ind, img in enumerate(img_tup,1) if any(img[sli])]
        nrows, ncols = img_tup[0][sli].shape
        base_item = np.zeros((nrows, ncols))
        for img, (start, stop) in zip(img_present, make_zebra_stripe_inds(0, ncols, len(img_present))):
            base_item[:,start:stop] = img[:,start:stop]
        baseimage[sli] = base_item
    
    return baseimage
        
    

# <codecell>

import colorsys

def get_colors(num_colors):
    colors=[]
    np.random.seed(50)
    for i in np.arange(0., 360., 360. / num_colors):
        hue = i/360.
        lightness = (50 + np.random.rand() * 10)/100.
        saturation = (90 + np.random.rand() * 10)/100.
        colors.append([255*num for num in colorsys.hls_to_rgb(hue, lightness, saturation)])
    return colors
    
def make_final_image(base_img, combined_image, num_images, colors = None):
    
    
    nums = set(combined_image.flat)
    nums.discard(0)
    if colors == None:
        colors = dict(zip(nums, get_colors(len(nums))))
    
    for num in nums:
        positions = np.argwhere(combined_image==num)
        for i,j in positions:
            base_img[i,j,:] = colors[num]
    return base_img
       
    

# <codecell>

import csv
from collections import defaultdict

def process_gene_list(gene_file, pathway_list):
    
    gene_dict = defaultdict(set)
    with open(gene_file) as handle:
        for row in csv.reader(handle, delimiter = '\t'):
            gene_dict[row[1]].add(row[0])
    groups = sorted(gene_dict.keys())
    for pathway in pathway_list:
        labeled_arrays = []
        baseimg = get_numpy_from_url(get_base_image_url(pathway))
        for key in groups:
            genes = gene_dict[key]
            nurl = get_img_url_from_genes(pathway, genes)
            narray = get_numpy_from_url(nurl)
            labeled, _ = get_colored_regions(narray)
            labeled_arrays.append(labeled)
        if len(labeled_arrays) == 0:
            continue
        combined = combine_images(labeled_arrays)
        try:
            nbase = make_final_image(baseimg.copy(), combined, len(labeled_arrays))
        except ZeroDivisionError:
            continue
        yield pathway, nbase, groups
        
    

# <codecell>

def make_legend(set_list, ax):
    
    nitems = len(set_list)
    colors = get_colors(nitems)
    ax.hold(True)
    for color, key in zip(colors, set_list):
        ax.plot([0],[0], label = key, color = np.array(color)/255)
    ax.legend()
    ax.hold(False)
    

# <codecell>

fname = '/home/will/HIVSystemsBio/out_gene_list.tsv'
plist = set(['hsa04060', 'hsa04062', 'hsa04064', 'hsa04620', 'hsa04623', 'hsa05132'])
for pathway, nbase, groups in process_gene_list(fname, sorted(plist)):
    print pathway
    figure(figsize=(20,20))
    ax = gca()
    ax.imshow(nbase)
    make_legend(groups, ax)
    plt.savefig('/home/will/HIVSystemsBio/figures/%s.png' % pathway)

# <codecell>


# <codecell>


# <codecell>


