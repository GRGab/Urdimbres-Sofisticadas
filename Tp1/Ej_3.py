#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  9 21:45:45 2018

@author: tomas
"""

import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from networkx.readwrite.gml import read_gml
from histograma import histograma
from __future__ import division

internet = read_gml('Tp1/tc01_data/as-22july06.gml')
nodes = []
degrees = []
for a, b in internet.degree():
    nodes.append(a)
    degrees.append(b)
#%% Comparación de visualizaciones

# Para comparar los bineados logarítmicos y no logarítmicos, lo justo es
# excluir a los nodos de grado 0 en ambos 

fig, axes = plt.subplots(4, 2, figsize=(10,10))
axes = axes.flatten()
logbinss = [0, 0, 0, 0, 1, 1, 1, 1]
logxs    = [0, 0, 1, 1, 0, 0, 1, 1]
logys    = [0, 1, 0, 1, 0, 1, 0, 1]

t = ['Bines lineales', 'Bines logarítmicos']
titulos  = [t[i] for i in logbinss]
xlabels = [('Grado (adim.)' if i in [6,7] else None) for i in range(8)]
ylabels = [(True if i % 2 == 0 else False) for i in range(8)]

for i in range(8):
    histograma(degrees,
               logbins=logbinss[i], logx=logxs[i], logy=logys[i], ax=axes[i],
               titulo=titulos[i], xlabel=xlabels[i], ylabel=ylabels[i],
               ecolor='k', errbars=False, 
               labelsize=10, ticksize=10,
               bins=(1, max(degrees) + 2, 100))
#%% Código viejo
    
#plt.hist([i for i in degree if i < 20], bins=20)
#plt.hist([i for i in degree if i > 20], bins=50)

k, values = np.histogram(degrees, bins = max(degrees))

p_k = [i/len(nodes) for i in k]

'''Algo interesante de ver es que el 98% de los degrees estan entre los 
degrees 0 y 20:    np.sum(p_k[0:20])'''
plt.figure(1)
plt.plot(p_k, '.')
plt.yscale('log')
plt.xscale('log')

plt.figure(2)
plt.hist(p_k, bins = 50)
plt.xscale('log')
plt.yscale('log')