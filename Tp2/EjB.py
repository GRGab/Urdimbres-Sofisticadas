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
from agregar_esencialidad import agregar_esencialidad

#%%
apms = ldata('Tp2/tc02Data/yeast_AP-MS.txt')

y2h = ldata('Tp2/tc02Data/yeast_Y2H.txt')

lit = ldata('Tp2/tc02Data/yeast_LIT.txt')

lit_r = ldata('Tp2/tc02Data/yeast_LIT_Reguly.txt')
lit_r = [fila[:2] for fila in lit_r]

g_apms = nx.Graph()
g_apms.add_edges_from(apms)

g_lit = nx.Graph()
g_lit.add_edges_from(lit)

g_lit_reg = nx.Graph()
g_lit_reg.add_edges_from(lit_r)

g_y2h = nx.Graph()
g_y2h.add_edges_from(y2h)

ess = ldata('Tp2/tc02Data/Essential_ORFs_paperHe.txt')
ess = [fila for fila in ess if (fila != [] and fila[0] != 'updated')]
ess = [fila[1] for fila in ess]

def k_medio(G):
    N = G.order()
    k_med = sum(k for (nodo, k) in G.degree) / N
    return k_med

#%%
################ Punto B:

######## Tabla 1 de Zotenko: propiedades de las bases de datos

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
######## Tabla 2 de Zotenko: overlap entre bases de datos

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
##### Figura 1.a de Zotenko

def ess_vs_k(g, ess):
    nodos, values = agregar_esencialidad(g, ess)
    nodos = []
    k = []
    for a, b in g.degree():
        nodos.append(a)
        k.append(b)
    k_norm = [i/max(k) for i in k]
    threshold = np.arange(0.01, 1, 0.001)
    y = []
    x = []
    for j in threshold:
        y_i = 0
        x_i = 0
        for i in range(len(nodos)):
            if k_norm[i]>j:
                if values[i] == 1:
                    y_i += 1
                x_i += 1
        if x_i != 0:
            y.append(y_i/x_i)
            x.append(x_i/len(nodos))
    return y, x


y_0, x_0 = ess_vs_k(g_apms, ess)
y_1, x_1 = ess_vs_k(g_lit, ess)
y_2, x_2 = ess_vs_k(g_y2h, ess)

plt.plot(x_0, y_0, label = 'Ap')
plt.plot(x_1, y_1, label = 'Lit')
plt.plot(x_2, y_2, label = 'Y2H')
plt.legend()




