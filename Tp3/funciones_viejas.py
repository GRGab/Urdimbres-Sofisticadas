# -*- coding: utf-8 -*-
"""
Created on Wed Oct 24 19:03:30 2018

@author: Gabo
"""

import matplotlib.pyplot as plt
import numpy as np

import time

import networkx as nx
from networkx.algorithms.community.centrality import girvan_newman
from networkx import NetworkXError
from networkx.utils import not_implemented_for
from networkx.algorithms.community.community_utils import is_partition

import igraph as igraph

import community # instalar como python-Louvain

def calcular_particion(np_adj_list, method):
    t0=time.time()
    if method in ['infomap', 'fastgreedy']:
        g = igraph.Graph.Weighted_Adjacency(np_adj_list.tolist(),mode=igraph.ADJ_UPPER)
        if method=="infomap":
            labels = g.community_infomap(edge_weights="weight").membership
    #    labels = g.community_label_propagation(weights="weight").membership
        if method=="fastgreedy":
            labels = g.community_fastgreedy(edge_weights="weight").membership
   if method == 'edge_bet':
       labels = girvan_newman(g)
    print("Duración: {}s".format(time.time()-t0))
    return labels

""" fast-greedy sin igraph"""
import networkx as nx
import matplotlib.pyplot as plt
G = nx.balanced_tree(h=3,r=2)
nx.draw(G,with_labels=True)
plt.show()

#%%

""" Hacer particion con louvain"""
partition = community.best_partition(G)
#%%
comus = nx.algorithms.community.greedy_modularity_communities(G, weight=None)
list(comus)
#%%
"""Girvan-newman"""
com = nx.algorithms.community.centrality.girvan_newman(G)
com = list(com)
com = [list([list(conjunto) for conjunto in tupla]) for tupla in com]
