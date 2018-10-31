#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 31 14:47:21 2018

@author: matias
"""

from percolation_clique import *
import networkx as nx
from matplotlib import pyplot as plt
from networkx.readwrite.gml import read_gml
from collections import defaultdict
import sys
sys.path.append('./Tp3/')
#%%
dolph = read_gml('Tp3/dolphins.gml')    

#genders = dict(ldata('Tp3/dolphinsGender.txt'))
#
## Agrego los sexos a los dicts de cada delfín
#for nodo, dict_nodo in dict(dolph.nodes).items():
#    dict_nodo['gender'] = genders[nodo]

colors = []

for j in [2, 3, 4, 5, 6, 7, 8]:
    c = list(k_clique_communities(dolph, j))
    for i in range(len(c)):
        c[i] = list(c[i])
        
    colores = comunidad_a_color(dolph, c)
    #le pongo color aqua a los que no tienen comunidad.
    for i in range(len(colores[0])):
        if colores[0][i] == 0.0:
            colores[0][i] = 'aqua'
            
    colors.append(colores[0])



#[ax1, ax2, ax3, ax4], [ax5, ax6, ax7, ax8]

ns = 35
#nx.draw(dolph, node_color = colors[0])
pos = nx.spring_layout(dolph)

fig, ([ax1, ax2, ax3], [ax4, ax5, ax6]) = plt.subplots(2, 3)

plt.sca(ax1)
ax1.set_title('k = 2')
nx.draw(dolph, node_size = ns, node_color=colors[0], pos = pos)

plt.sca(ax2)
ax2.set_title('k = 3')
nx.draw(dolph, node_size = ns, node_color=colors[1], pos = pos)

plt.sca(ax3)
ax3.set_title('k = 4')
nx.draw(dolph, node_size = ns, node_color=colors[2], pos = pos)

plt.sca(ax4)
ax4.set_title('k = 5')
nx.draw(dolph, node_size = ns, node_color=colors[3], pos = pos)

plt.sca(ax5)
ax5.set_title('k = 6')
nx.draw(dolph, node_size = ns, node_color=colors[4], pos = pos)

plt.sca(ax6)
ax6.set_title('k = 7')
nx.draw(dolph, node_size = ns, node_color=colors[5], pos = pos)
