# -*- coding: utf-8 -*-
"""
Created on Fri Oct  5 18:00:36 2018

@author: Gabo
"""
import networkx as nx
import numpy as np

#Cambio de Mati
def agregar_esencialidad(G, ess):
    g = list(G.nodes())
    value = np.zeros([len(g)])
    nodos_esenciales = []
    for h in range(len(g)):
        for i in range(2,1158):
            if g[h] == ess[i][1]:
                value[h] = 1
                nodos_esenciales.append(g[h])
                break
    return nodos_esenciales, value

#%%
            
#Funciones de tomi
#import networkx as nx
#import numpy as np
#def agregar_esencialidad(G, ess):
#    g = list(G.nodes())
#    value = np.zeros([len(g)])
#    for h in range(len(g)):
#        for i in range(len(ess)):
#                if g[h] == ess[i]:
#                    value[h] = 1
#                    break
#    return g, value
#
#def agregar_esencialidad_dict(G, esenciales):
#    """ Modifica inplace el grafo G"""
#    for nodo, dict_nodo in dict(G.nodes).items():
#        if nodo in esenciales:
#            dict_nodo['esencialidad'] = 1
#        else:
#            dict_nodo['esencialidad'] = 0