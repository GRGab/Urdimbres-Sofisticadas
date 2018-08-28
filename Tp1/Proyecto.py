#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 24 18:06:47 2018

@author: matias
"""

from lectura import ldata
import networkx as nx
from matplotlib import pyplot as plt

apms = ldata('/home/matias/Documentos/Facultad/Redes/Urdimbres-Sofisticadas/'
             +'Tp1/tc01_data/yeast_AP-MS.txt')

lit = ldata('/home/matias/Documentos/Facultad/Redes/Urdimbres-Sofisticadas/'
             +'Tp1/tc01_data/yeast_LIT.txt')

y2h = ldata('/home/matias/Documentos/Facultad/Redes/Urdimbres-Sofisticadas/'
             +'Tp1/tc01_data/yeast_Y2H.txt')

#%%
def es_dirigido(data):
    """data debe ser una lista de enlaces (tuplas). Si el resultado es 0,
    entonces es no dirigido; si es distinto de cero, es dirigido."""
    n = 0
    for (x, y) in data:
        for (a, b) in data:
            if a == y and b == x:
                n += 1
    return n/2
print(es_dirigido(apms), es_dirigido(lit),
      es_dirigido(y2h))
# Vemos que y2h y lit son dirigidos, mientras que apms no.
#%%

g_apms = nx.Graph()
g_apms.add_edges_from(apms)
#nx.draw(g_apms, node_size = 50)
#
g_lit = nx.DiGraph()
g_lit.add_edges_from(lit)
#nx.draw(g_lit, node_size = 50)
#
g_y2h = nx.DiGraph()
g_y2h.add_edges_from(y2h)
#nx.draw(g_y2h, node_size = 35)
#%%
f, (ax1, ax2, ax3) = plt.subplots(1, 3)
plt.sca(ax1)
ax1.set_title('Lit')
nx.draw(g_lit, node_size = 10)

plt.sca(ax2)
ax2.set_title('Y2h')
nx.draw(g_y2h, node_size = 10)

plt.sca(ax3)
ax3.set_title('APMS')
nx.draw(g_apms, node_size = 10)

#%%
#Cuantos nodos hay?
print('El número de nodos de cada grafo es',
      g_lit.order(), g_y2h.order(), g_apms.order()) 
#Enlaces
print('El número de enlaces para cada grafo es',
      g_lit.size(), g_y2h.size(), g_apms.size())
#%%
def k_medio(G):
    N = G.order()
    if isinstance(G, nx.DiGraph):
        kin_med = sum(k for (nodo, k) in G.in_degree) / N
        kout_med = sum(k for (nodo, k) in G.out_degree) / N
        return kin_med, kout_med
    else:
        k_med = sum(k for (nodo, k) in G.degree) / N
        return k_med
        
kin_medio_lit, kout_medio_lit = k_medio(g_lit)
kin_medio_y2h, kout_medio_y2h = k_medio(g_y2h)
k_medio_apms = k_medio(g_apms)
print(kin_medio_lit, kout_medio_lit,
      kin_medio_y2h, kout_medio_y2h,
      k_medio_apms)

#Valor de k max y min
def k_extremos(G):
   k_min = min(k for (nodo, k) in G.degree)
   k_max = max(k for (nodo, k) in G.degree)
   return k_min, k_max
k_min_apms, k_max_apms = k_extremos(g_apms)
k_min_lit, k_max_lit = k_extremos(g_lit)
k_min_y2h, k_max_y2h = k_extremos(g_y2h)
print(k_min_apms, k_max_apms)
print(k_min_y2h, k_max_y2h)
print(k_min_lit, k_max_lit)

#%%
#Coeficientes de clustering <C_i> y C_Δ de la red.
#C_Δ
print (nx.transitivity(g_apms),nx.transitivity(g_y2h),nx.transitivity(g_lit))
#<C_i> (No definido para apms y y2h. Mati: pensar porque es asi)
print(nx.average_clustering(g_apms))

#%%
#Densidad de la red
print('La desnidad de las redes son', nx.density(g_lit), nx.density(g_y2h)
, nx.density(g_apms))
#%%
#Diámetro de la red 
#(Mati: tira error-->'infinito path lenght porque el diagrama no esta
# fuertemente conectado.' Hay que consultar.)

print('El diámetro de las redes son', nx.diameter(g_lit,e=None)
,nx.diameter(g_y2h,e=None), nx.diameter(g_apms,e=None))










