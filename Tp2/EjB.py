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


import sys
sys.path.append('./Tp2/')
from agregar_esencialidad import agregar_esencialidad, agregar_esencialidad_dict

#%%
apms = ldata('Tp2/tc02Data/yeast_AP-MS.txt')

y2h = ldata('Tp2/tc02Data/yeast_Y2H.txt')

lit = ldata('Tp2/tc02Data/yeast_LIT.txt')

lit_r = ldata('Tp2/tc02Data/yeast_LIT_Reguly.txt')
lit_r = [fila[:2] for fila in lit_r[1:]]

g_apms = nx.Graph()
g_apms.add_edges_from(apms)

g_lit = nx.Graph()
g_lit.add_edges_from(lit)

g_lit_reg = nx.Graph()
g_lit_reg.add_edges_from(lit_r)

g_y2h = nx.Graph()
g_y2h.add_edges_from(y2h)

ess = ldata('Tp2/tc02Data/Essential_ORFs_paperHe.txt')
ess =  ess[2:-4]
ess = [fila[1] for fila in ess]
ess = np.unique(ess)

for g in [g_apms, g_lit, g_lit_reg, g_y2h]:
    agregar_esencialidad_dict(g, ess)

def k_medio(g):
    N = g.order()
    k_med = sum(k for (nodo, k) in g.degree) / N
    return k_med

#%%
################ Punto B:

######## Tabla 1 de Zotenko: propiedades de las bases de datos

def fila_tabla_1(g):
    g = max(nx.connected_component_subgraphs(g), key=len)
    N = g.order()
    L = g.size()
    K = k_medio(g)
    Ci = np.average(list(nx.clustering(g).values()))
    return N, L, K, Ci

grafos = [g_lit, g_lit_reg, g_apms, g_y2h]
Ns, Ls, Ks, Cis = [], [], [], []
for g in grafos:
    N, L, K, Ci = fila_tabla_1(g)
    Ns.append(N)
    Ls.append(L)
    Ks.append(K)
    Cis.append(Ci)

tabla1 = pd.DataFrame(data={'Número de nodos':Ns,
                            'Número de enlaces':Ls,
                            'Grado medio':Ks,
                            'Coef. de clustering medio':Cis}, 
                      index=['Lit', 'Lit_reg', 'APMS', 'Y2H'])
pd.set_option('precision', 2)
tabla1 

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

listas_enlaces = [lit, lit_r, apms, y2h]
enlaces_compartidos = np.zeros((4,4))
for i in range(4):
    for j in range(4):
        if i==j:
            enlaces_compartidos[i,j] = 1
        else:
            enlaces_compartidos[i,j] = comparten_enlaces(listas_enlaces[i],
                                                         listas_enlaces[j])
            
print(enlaces_compartidos)
# Vemos que este método falla, asignando 373% de enlaces compartidos con
# lit_reg en la red lit.
# Esto probablemente tiene que ver con que este método no considera correctamente
# la posibilidad de existencia de repetidos en las listas de enlaces.
# Además, habría que hacer la cuenta únicamente sobre la componente gigante.
# Para esto nos conviene trabajar con los grafos directamente.
#%%   

def comparten_enlaces2(grafo1, grafo2):
    """Por convención, devuelve la fracción de enlaces del grafo 1 que se 
    encuentran en el grafo 2. """
    # Nos quedamos con las componentes gigantes
    grafo1 = max(nx.connected_component_subgraphs(grafo1), key=len)
    grafo2 = max(nx.connected_component_subgraphs(grafo2), key=len)
    # Comienza el cálculo
    n = 0
    for (x, y) in grafo1.edges:
        for (a, b) in grafo2.edges:
            if (a == y and b == x) or (a == x and b == y):
                n += 1
                break # Ya no va a encontrar más así que nos ahorramos pasos
    return n/grafo1.size()

grafos = [g_lit, g_lit_reg, g_apms, g_y2h]
enlaces_compartidos2 = np.zeros((4,4))
for i in range(4):
    for j in range(4):
        if i==j:
            enlaces_compartidos2[i,j] = 1
        else:
            enlaces_compartidos2[i,j] = comparten_enlaces2(grafos[i], grafos[j])

#%% Comparamos ambos métodos
np.savez('Tp2/tc02Data/tabla2zotenko.npz', data=enlaces_compartidos2)
print(enlaces_compartidos2)

# Este segundo método da resultados razonables en todos los casos.
#%%
tabla2 = pd.DataFrame()
#%%
##### Figura 1.a de Zotenko

def essfrac_vs_hubfrac(g):
    """Toma grafo que ya tiene incorporada la información de esencialidad
    Trabaja únicamente sobre la componente gigante (como hicieron Zotenko et al.)"""
    
    g = max(nx.connected_component_subgraphs(g), key=len)
    kmax = max(dict(g.degree).values())
    num_hubs, num_hubs_esens = np.zeros((2, kmax+1))
    for node, degree in g.degree():
        num_hubs[:degree+1] += np.ones((degree+1))
        if g.nodes[node]['esencialidad'] == 1:
            num_hubs_esens[:degree+1] += np.ones((degree+1))
    frac_hubs_sonesens = num_hubs_esens / num_hubs
    frac_sonhubs = num_hubs / g.order()
    return frac_sonhubs, frac_hubs_sonesens


#%%
        
x_0, y_0 = essfrac_vs_hubfrac(g_apms)
x_1, y_1 = essfrac_vs_hubfrac(g_lit)
x_2, y_2 = essfrac_vs_hubfrac(g_y2h)
x_3, y_3 = essfrac_vs_hubfrac(g_lit_reg)

#%%
fontsize = 18
ticksize = 16
with plt.style.context(('seaborn')):
    fig, ax = plt.subplots(figsize=(12,8))
ax.plot(x_1, y_1, '-o', label = 'Lit')
ax.plot(x_3, y_3, '-o', label = 'Lit_reg')
ax.plot(x_0, y_0, '-o', label = 'APMS')
ax.plot(x_2, y_2, '-o', label = 'Y2H')
ax.legend(fontsize=fontsize)
ax.tick_params(labelsize=ticksize)
ax.set_xlabel('Fracción de nodos considerados como hubs',
              fontsize=fontsize)
ax.set_ylabel('Fracción de hubs que son esenciales',
              fontsize=fontsize)
fig.tight_layout()

#%% Versión vieja  de la figura 1.a Zotenko (Tomi)

# =============================================================================
# def ess_vs_k(g, ess):
#     _, values = agregar_esencialidad(g, ess)
#     nodos = []
#     k = []
#     for a, b in g.degree():
#         nodos.append(a)
#         k.append(b)
#     k_norm = [i/max(k) for i in k]
#     threshold = np.arange(0.01, 1, 0.001)
#     y = []
#     x = []
#     for j in threshold:
#         y_i = 0 # Núm. de hubs que son esenciales
#         x_i = 0 # Núm. de nodos que son considerados hubs
#         for i in range(len(nodos)):
#             if k_norm[i]>j:
#                 if values[i] == 1:
#                     y_i += 1
#                 x_i += 1
#         if x_i != 0:
#             y.append(y_i/x_i)
#             x.append(x_i/len(nodos))
#     return y, x
# =============================================================================

# =============================================================================
# y_0, x_0 = ess_vs_k(g_apms, ess)
# y_1, x_1 = ess_vs_k(g_lit, ess)
# y_2, x_2 = ess_vs_k(g_y2h, ess)
# y_3, x_3 = ess_vs_k(g_lit_reg, ess)
# 
# plt.figure()
# plt.plot(x_0, y_0, label = 'Ap')
# plt.plot(x_1, y_1, label = 'Lit')
# plt.plot(x_2, y_2, label = 'Y2H')
# plt.plot(x_3, y_3, label = 'Lit_r')
# plt.legend()
# 
# =============================================================================
