# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 18:52:38 2018

@author: Gabo
"""

import networkx as nx
from networkx.algorithms.community.centrality import girvan_newman
import community # Python-Louvain
import matplotlib.pyplot as plt
#%%

G = nx.karate_club_graph()
plt.figure()
nx.draw(G,with_labels=True)
plt.show()

partition = community.best_partition(G)
partition