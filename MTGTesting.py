# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import os, os.path
from matplotlib import pyplot as plt
from pylab import get_cmap
import SimpleCV as cv
from glob import glob
from urllib import urlopen
from bs4 import BeautifulSoup

# <codecell>

data = urlopen('http://deckbox.org/games/mtg/cards?p=0')
soup = BeautifulSoup(data)

# <codecell>

table = soup.find(attrs = {'class':'tabular_cards'})
for row in table.find_all('tr'):
    cols = row.find_all('td')
    if len(cols) > 2:
        print cols[0].text.strip(), cols[-1].text.strip()

# <codecell>

cols[-1].text.strip()

# <codecell>


