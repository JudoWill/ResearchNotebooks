# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import networkx as nx
import matplotlib.pyplot as plt
from pandas import *
import csv
import os, os.path

os.chdir('/home/will/Dropbox/CircusInstruction/')

# <codecell>

obj = ExcelFile('TrickNetwork.xls')
trick_data = obj.parse('Trick Names', index_col = 0)
trick_links = obj.parse('Trick Links')

# <codecell>

G = nx.DiGraph()
G.add_nodes_from(trick_data.index)

# <codecell>

edge_color_types = {
'Drop':'R',
'Roll':'C',
'Wrap':'B',
'Unwrap':'B',
'Twist':'C',
'Invert':'M',
'Side Lift':'Y',
'Climb Over':'G',
'Swing':'k'
}

edge_colors = {}


for num, row in trick_links.iterrows():
    if (row['Source'] not in trick_data.index) or (row['Destination'] not in trick_data.index):
        print 'Wrong in line', num+1, row
    else:
        G.add_edge(row['Source'], row['Destination'])
        edge_colors[(row['Source'], row['Destination'])] = edge_color_types.get(row['Transition'], 'k')
        if row['Reversible']:
            G.add_edge(row['Destination'], row['Source'])

# <codecell>



typ_colors = {
'Drop':'R',
'Climb':'G',
'Wrap':'B',
'Pose':'Y',
'Transition':'C'
}
node_colors = {}
for name, typ in zip(trick_data.index, trick_data['Type'].values):
    node_colors[name] = typ_colors[typ]
    
node_sizes = nx.degree(G)


#plt.savefig('TrickNetwork.png')

# <codecell>



plt.figure(figsize = (40,40))
nx.draw_graphviz(G, alpha = 0.2, 
        node_color = [node_colors[n] for n in G.nodes()],
        node_size = [500*node_sizes[n] for n in G.nodes()],
        edge_color = [edge_colors.get(e, 'k') for e in G.edges()],
        scale = 500, iterations = 5000)

plt.savefig('TrickNetwork-tmp.png')
plt.close()

# <codecell>

trap_network = G.subgraph(trick_data[trick_data['Trapeze']==1].index)
node_sizes = nx.degree(trap_network)
plt.figure(figsize = (40,40))
nx.draw_graphviz(trap_network, alpha = 0.2, 
        node_color = [node_colors[n] for n in trap_network.nodes()],
        node_size = [500*node_sizes[n] for n in trap_network.nodes()],
        edge_color = [edge_colors.get(e, 'k') for e in trap_network.edges()],
        scale = 500, iterations = 5000)

plt.savefig('Trap-TrickNetwork.png')
plt.close()

# <codecell>

sum(1/x for x in dists.values() if x)/len(dists)

# <codecell>

dists

# <codecell>


