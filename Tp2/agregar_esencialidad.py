# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 18:00:36 2018

@author: Gabo
"""
import networkx as nx
import numpy as np

def agregar_esencialidad(G, ess):
    g = list(G.nodes())
    value = np.zeros([len(g)])
    for h in range(len(g)):
        for i in range(len(ess)):
            for j in range(len(ess[i])):
                if g[h] == ess[i][j]:
                    value[h] = 1
                    break
    return g, value