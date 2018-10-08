# -*- coding: utf-8 -*-
            
import networkx as nx
import numpy as np
def agregar_esencialidad(G, ess):
    g = list(G.nodes())
    value = np.zeros([len(g)])
    for h in range(len(g)):
        for i in range(len(ess)):
                if g[h] == ess[i]:
                    value[h] = 1
                    break
    return g, value

def agregar_esencialidad_dict(G, esenciales):
    """ Modifica inplace el grafo G"""
    for nodo, dict_nodo in dict(G.nodes).items():
        if nodo in esenciales:
            dict_nodo['esencialidad'] = 1
        else:
            dict_nodo['esencialidad'] = 0