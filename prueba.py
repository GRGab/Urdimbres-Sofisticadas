# -*- coding: utf-8 -*-
"""
Created on Fri Aug 24 18:33:13 2018

@author: Gabo
"""

import networkx as nx
import matplotlib.pyplot as plt

g = nx.Graph()
g.add_edges_from([(1,2), (1,3), (1,4), (3,4)])
h = nx.Graph()
h.add_edges_from([(1,2), (1,3), (1,4), (4,5), (5,1)])

fig, (ax1, ax2) = plt.subplots(1, 2)
plt.sca(ax1)
plt.title("Una red")
nx.draw(g)
plt.sca(ax2)
plt.title("Otra red")
nx.draw(h)
