# -*- coding: utf-8 -*-
"""
Created on Fri Aug 24 18:33:13 2018

@author: Gabo
"""

import networkx as nx
import matplotlib.pyplot as plt
import sys
sys.path.append('./Tp1/')
from dolphins_funciones import (contar_enlaces_internos,
                                contar_enlaces_entre_grupos)

g = nx.Graph()
g.add_edges_from([(1,2), (1,3), (1,4), (3,4)])
h = nx.Graph()
h.add_edges_from([(1,2), (1,3), (1,4), (4,5), (5,1)])

fig, (ax1, ax2) = plt.subplots(1, 2)
plt.sca(ax1)
plt.title("Una red")
nx.draw(g, with_labels=True)
plt.sca(ax2)
plt.title("Otra red")
nx.draw(h, with_labels=True)

#%% Clustering

clust_g = nx.clustering(g)
print(clust_g)
clust_h = nx.clustering(h)
print(clust_h)


#%% Contar enlaces internos

for nodo, dictattr in dict(h.nodes()).items():
    dictattr['aroma'] = 'rico' if nodo in [1,4,5] else 'feo'

#print(contar_enlaces_internos(h, 'aroma', 'rico'))
print(contar_enlaces_entre_grupos(h, 'aroma'))