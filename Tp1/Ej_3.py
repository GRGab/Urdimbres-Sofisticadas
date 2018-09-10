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
from lectura import ldata

internet = read_gml('/home/tomas/Desktop/Redes complejas/Urdimbres-Sofisticadas/Tp1/tc01_data/as-22july06.gml')
nodes = []
degree = []
for a, b in internet.degree():
    nodes.append(a)
    degree.append(b)
    
#nx.draw(internet)

#plt.hist([i for i in degree if i < 20], bins=20)
#plt.hist([i for i in degree if i > 20], bins=50)

k, values = np.histogram(degree, bins = max(degree))

p_k = [i/float(len(nodes)) for i in k]

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