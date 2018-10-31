#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 29 22:26:19 2018

@author: matias
"""

import networkx as nx
import numpy as np
from matplotlib import pyplot as plt
import matplotlib
import pandas as pd
import igraph as igraph
import time
from lectura import ldata

from os.path import join as osjoin
import sys
sys.path.append('./Final/')
from histograma import histograma
sys.path.append('./Tp3/')
from funciones_tp3 import calcular_particion, comunidad_a_color

#from funciones_ej_a(Mati) import ldata, calculate_infomap, modularity
#
#from scipy.cluster.hierarchy import dendrogram, linkage, cophenet
#from scipy.spatial.distance import pdist


def comunidad_a_color(g, lista):
    """
    Funcion para asignar colores a las comunidades. Devuelve una lista con colores
    con el mismo orden que la lista de nodos, para meter en la funcion 
    nx.draw(G, node_color = comunidad_a_color(G, lista)). 
    La lista corresponde a la lista de listas, donde cada sublista corresponde a
    una comunidad.
    
    Input: (g, lista)   (grafo de networkx, lista de listas)
    
    Returns:
            colores   (lista)
    .
    .
    """
    colores_posibles =  [b for (a,b) in matplotlib.colors.cnames.items()]
    colores_random = np.random.choice(np.arange(len(colores_posibles)), size=len(lista),
                                      replace=False)
    nodos = list(g.nodes())
    colores = list(np.zeros(len(nodos)))
    for i in range(len(lista)):
        for j in range(len(nodos)):
            index = colores_random[i]
            if nodos[j] in lista[i]:
                colores[j] = colores_posibles[index]
    return colores


#%% Punto 1    
wiki = ldata('Final/links.tsv')
g_wiki = nx.Graph()
g_wiki.add_edges_from(wiki)
#%%
wiki_infomap = calcular_particion(g_wiki)
wiki_infomap_color = comunidad_a_color(g_wiki, wiki_infomap)
wiki_fastgreedy = calcular_particion(g_wiki, method = 'fastgreedy')
wiki_fastgreedy_color = comunidad_a_color(g_wiki, wiki_fastgreedy)
#%%Ploteamos la red
fig, ax = plt.subplots()
plt.sca(ax)
ax.set_title('Wikipedia')
#nx.draw(g_wiki, node_color = wiki_infomap_color, node_size = 10)
nx.draw_spring(g_wiki, node_color = wiki_fastgreedy_color,
               node_size = 10, with_labels=True)
