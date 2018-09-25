#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 24 18:06:47 2018

@author: matias
"""

from lectura import ldata
import networkx as nx
import numpy as np
from matplotlib import pyplot as plt


apms = ldata('Tp1/tc01_data/yeast_AP-MS.txt')

lit = ldata('Tp1/tc01_data/yeast_LIT.txt')

y2h = ldata('Tp1/tc01_data/yeast_Y2H.txt')

#%%
def es_dirigido(data):
    """data debe ser una lista de enlaces (tuplas). Si el resultado es 0,
    entonces es no dirigido; si es distinto de cero, es dirigido.
    
    Esta función implementa una manera de inferir si un grafo es dirigido o no
    a partir de cómo se presentan sus enlaces en una lista; sin embargo este
    criterio no es infalible. De hecho, nos fue revelado que las 3 redes
    consideradas en este ejercicio son efectivamente no dirigidas, lo cual
    tiene sentido pues se trata de interacciones entre proteínas en las cuales
    no es evidente qué significado podría tener la direccionalidad."""
    n = 0
    for (x, y) in data:
        for (a, b) in data:
            if a == y and b == x:
                n += 1
    return n/2

#%%Iportamos y graficamos por separado cada red
g_apms = nx.Graph()
g_apms.add_edges_from(apms)
#plt.figure()
#nx.draw(g_apms, node_size = 50)
#
#g_lit = nx.DiGraph()
g_lit = nx.Graph()
g_lit.add_edges_from(lit)
#plt.figure()
#nx.draw(g_lit, node_size = 50)
#
#g_y2h = nx.DiGraph()
g_y2h = nx.Graph()
g_y2h.add_edges_from(y2h)
#plt.figure()
#nx.draw(g_y2h, node_size = 35)
#%% Graficamos en subplots
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
#%% Son dirigidas?
print(es_dirigido(apms), es_dirigido(lit),
      es_dirigido(y2h))
# Vemos que y2h y lit son dirigidos, mientras que apms no.
#%%
print(es_dirigido(apms) / g_apms.size(), es_dirigido(lit) / g_lit.size(),
      es_dirigido(y2h) / g_y2h.size())

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
    else:
        kin_med, kout_med = 0, 0
    k_med = sum(k for (nodo, k) in G.degree) / N
    return kin_med, kout_med, k_med
        
kin_medio_lit, kout_medio_lit, k_medio_lit = k_medio(g_lit)
kin_medio_y2h, kout_medio_y2h, k_medio_y2h = k_medio(g_y2h)
_, _, k_medio_apms = k_medio(g_apms)

print('Grados medios')
print(kin_medio_lit, kout_medio_lit, k_medio_lit)
print(kin_medio_y2h, kout_medio_y2h, k_medio_y2h)
print(k_medio_apms)

# kin y kout medios dan lo mismo. Parece no trivial pero es trivial ! :D ;)

#Valor de k max y min
def k_extremos(G):
   k_min = min(k for (nodo, k) in G.degree)
   k_max = max(k for (nodo, k) in G.degree)
   return k_min, k_max
k_min_apms, k_max_apms = k_extremos(g_apms)
k_min_lit, k_max_lit = k_extremos(g_lit)
k_min_y2h, k_max_y2h = k_extremos(g_y2h)
print('Grados extremos')
print(k_min_apms, k_max_apms)
print(k_min_y2h, k_max_y2h)
print(k_min_lit, k_max_lit)

#%%
#Densidad de la red
print('La densidad de las redes es', nx.density(g_lit), nx.density(g_y2h),
      nx.density(g_apms))
#%%
#Coeficientes de clustering <C_i> y C_Δ de la red.
print('C_Δ')
print (nx.transitivity(g_apms), nx.transitivity(g_y2h),
       nx.transitivity(g_lit))
# Ignora diferencia entre in y out
#<C_i> (No definido para apms y y2h. Mati: pensar porque es asi)
print('<C_i>')
def clustering(nodo):
    # Esto sería para grafos dirigidos
    pass
def clustering_medio(G):
    # Esto sería para grafos dirigidos
    return np.average(list(dict(nx.clustering(G)).values()))

print(clustering_medio(g_apms))
print(clustering_medio(g_lit))

#%%
#Diámetro de la red 

print('El diámetro de las redes es', nx.diameter(g_lit,e=None),
nx.diameter(g_y2h,e=None), nx.diameter(g_apms,e=None))
#%% Separamos la componente gigante de la red
nx.connected_component_subgraphs(g_apms, copy=True)

#%%
# Pregunta extra sobre APMS
# Pareciera estar hecho con el método de conectar
# todas las proteínas con todas porque tiene densidad y coefs de clustering
# mayores que los otros dos datasets, los cuales a su vez tienen valores
# relativamente similares para estas magnitudes

# Por otro lado, se observa mirando a ojo las componentes más pequeñas de la
# red que no parecen ser cliques, sino que tienen nodos especiales que funcionan
# de hubs, además de tener subgrafos que sí parecen ser cliques.
# Incluso en la componente gigante, parece muy preponderante esta estructura
# "arbórea" en la que muchos nodos "salen" de un único nodo, y después se
# conectan más o menos densamente entre sí.

# Estas observaciones parecerían refutar que la red esté hecha con el método de conectar
# todas con todas.