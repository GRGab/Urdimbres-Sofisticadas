import matplotlib.pyplot as plt
import numpy as np

import time

import networkx as nx
from networkx.algorithms.community.centrality import girvan_newman
from networkx import NetworkXError
from networkx.utils import not_implemented_for
from networkx.algorithms.community.community_utils import is_partition
from networkx.readwrite.gml import read_gml
from lectura import ldata

import sys
sys.path.append('./Tp3/')

import igraph as igraph
import community # instalar como python-Louvain

#%%
dolph = read_gml('dolphins.gml')
genders = dict(ldata('dolphinsGender.txt'))

def comunidad_a_color(g, lista):
    """
    Funcion para asignar colores a las comunidades. Devuelve una lista con colores
    con el mismo orden que la lista de nodos, para meter en la funcion 
    nx.draw(G, node_color = comunidad_a_color(G, lista)). 
    La lista corresponde a la lista de listas, donde cada sublista corresponde a
    una comunidad.
    
    Input: (g, lista)   (grafo de networkx, lista de listas)
    
    Returns:
            colores   (lista)
    .
    .
    """
    colores_posibles = ['r', 'b', 'g', 'k', 'c', 'y', 'violet', 'sandybrown', 'orange', 'indianred', 'darkgray', 'darksalmon']
    colores_random = np.random.randint(len(colores_posibles), size = len(lista))
    nodos = list(g.nodes())
    colores = list(np.zeros(len(nodos)))
    for i in range(len(lista)):
        for j in range(len(nodos)):
            index = colores_random[i]
            if nodos[j] in lista[i]:
                colores[j] = colores_posibles[index]
    return colores
