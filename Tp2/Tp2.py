#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 28 18:34:25 2018

@author: matias
"""

from lectura import ldata
import networkx as nx
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd

#%%
apms = ldata('tc02Data/yeast_AP-MS.txt')

y2h = ldata('tc02Data/yeast_Y2H.txt')

lit = ldata('tc02Data/yeast_LIT.txt')

lit_r = ldata('tc02Data/yeast_LIT_Reguly.txt')


g_apms = nx.Graph()
g_apms.add_edges_from(apms)

g_lit = nx.Graph()
g_lit.add_edges_from(lit)
#
#g_lit_reg = nx.Graph()
#g_lit_reg.add_edges_from(lit_r)

g_y2h = nx.Graph()
g_y2h.add_edges_from(y2h)

ess = ldata('tc02Data/Essential_ORFs_paperHe.txt')




#%%

def k_medio(G):
    N = G.order()
    k_med = sum(k for (nodo, k) in G.degree) / N
    return k_med



def tabla_1(red):
    red = max(nx.connected_component_subgraphs(red), key=len)
    N = red.order()
    L = red.size()
    K = k_medio(red)
    Ci = np.average(list(nx.clustering(red).values()))
    return [N,L,K,Ci]

data = pd.DataFrame({"Nombre de la red": ['Y2H','AP-MS','Lit'],
                     "N":[tabla_1(g_y2h)[0],tabla_1(g_lit)[0],tabla_1(g_apms)[0]],
                     "L":[tabla_1(g_y2h)[1],tabla_1(g_lit)[1],tabla_1(g_apms)[1]],
                     "$$<C_{i}>$$":[tabla_1(g_y2h)[2],tabla_1(g_lit)[2],tabla_1(g_apms)[2]],
                    })#empty dataframe
data
#%%
#Punto b
#Comparamos las listas de enlaces

def comparten_enlaces(data1,data2):
    """data1 y data2 son listas de enlaces (tuplas) que voy a comparar. Por 
    convención, voy a comparar data1 con data2 (no es conmutativo!). El
    resultado en la salida es la proporcion de enlaces de data 1 que se 
    encuentran en data2"""
    n = 0
    for (x, y) in data1:
        for (a, b) in data2:
            if (a == y and b == x) or (a == x and b == y):
                n += 1
    g = nx.Graph()
    g.add_edges_from(data1)
    return n/g.size()

comparten_enlaces(apms,apms)
#%%
def lectura(archive):
    f = open(archive)
    data = []
    for line in f:
        line = line.strip()
        col = line.split()
#        import pdb; pdb.set_trace()
        data.append(col)	
    return data

ess = lectura('tc02Data/Essential_ORFs_paperHe.txt')








