#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 28 18:37:39 2018

@author: matias
"""

import numpy as np
from networkx.readwrite.gml import read_gml
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import networkx as nx
from networkx.algorithms.community.community_utils import is_partition

from networkx.readwrite.gml import read_gml
import sys
sys.path.append('./Tp3')
from funciones_tp3 import (calcular_particion, NotAPartition, indices_to_nodos_particion,
                           comunidad_a_color, graficar_dist)
from guardar_particiones import guardar_particiones
#%%
dolph = read_gml('Tp3/dolphins.gml')    
lista_de_metodos = ["infomap","label_prop", "fastgreedy", "eigenvector",
                    "louvain", "edge_betweenness", "walktrap"]
guardar_particiones(dolph, 200,1000, lista_de_metodos)
#%% Importamos todo y graficamos
npzfile = np.load('Tp3/tc03Data/Ej_b_particiones_numeros.npz')
rewire = npzfile['mod_rewire']
original = npzfile['mod_original']
for i in range(0,7):
    graficar_dist(dolph, lista_de_metodos, rewire[i], original[i],
                  metodo = i)