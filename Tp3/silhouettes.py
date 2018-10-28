# -*- coding: utf-8 -*-
"""
Created on Wed Oct 24 19:38:39 2018

@author: Gabo
"""

import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from networkx.algorithms.community.community_utils import is_partition

from networkx.readwrite.gml import read_gml


import sys
sys.path.append('./Tp3')
from funciones_tp3 import calcular_particion, NotAPartition, indices_to_nodos_particion

from networkx.readwrite.gml import read_gml
#%%
def donde(nodo, particion):
    """Devuelve el índice del cluster de la partición al que pertence el nodo."""
    for j, cluster in enumerate(particion):
        if nodo in cluster:
            return j

def silhouettes(G, particion):
    """
    Calcula el valor de silhouette para cada nodo del grafo 'G' dada una
    partición 'particion' como lista de listas. Dicho valor está dado por
    
    s(i) = (b(i) - a(i)) / max(a(i), b(i))
    
    donde a(i) es la distancia media a todos los nodos del mismo cluster que i
    y b(i) es la mínima de las distancias medias a los distintos clusters a los
    cuales no pertenece i. Para mayor claridad, sea c_i el cluster al que
    pertenece i, y sea Q = particion - c_i el conjunto de los clusters a los cuales
    no pertenece i. Entonces se define
    
    b(i) = min{promedio{d(i,j) : j in cluster} : cluster in Q}
    
    b(i) también se suele llamar "distancia media al cluster más cercano".
    
    'particion' es lista de listas. Cada sublista es un cluster y sus elementos
    son los nombres de los nodos que pertenecen a dicho cluster.
    """
    
    if not is_partition(G, particion):
        raise NotAPartition(G, particion)

    ds = list(nx.all_pairs_shortest_path_length(G))
    d = lambda i,j: ds[i][1][j]
    # ds[i][1][j] es la distancia (longitud del camino más corto)
    # entre i y j
    
    n = G.order()
    output = np.zeros((n))
    for i, nodo in enumerate(G.nodes()):
        k = donde(nodo, particion)
        cluster_actual = particion[k]
        otros_clusters = (particion[j] for j in range(len(particion)) if j != k)
        a = np.average([d(i,j) for j in cluster_actual])
        
        dists_interclusters = [np.average([d(i,j) for j in cluster]) \
                               for cluster in otros_clusters]
        b = min(dists_interclusters)
        
        output[i] = b - a / max(a, b)
    return output

#def graficar_silhouettes():
    
    
#%%
if __name__ == '__main__':
    G = nx.balanced_tree(h=3,r=2)
    particion = calcular_particion(G, method='infomap')
    sil = silhouettes(G, particion)
    #%%
    dolph = read_gml('Tp3/dolphins.gml')
    plt.figure(); nx.draw(dolph)
    
    npzfile = np.load('Tp3/tc03Data/Ej_b_particiones.npz')
    rewire = npzfile['salida']
    original = npzfile['salida_grafo_original']
    part = original[0]
    sil = silhouettes(dolph, part)
        