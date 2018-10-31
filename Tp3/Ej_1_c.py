#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 31 13:55:38 2018

@author: matias
"""

import numpy as np
from networkx.readwrite.gml import read_gml
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import networkx as nx
import pandas as pd
from networkx.algorithms.community.community_utils import is_partition
import sys
sys.path.append('./Tp3')
from lectura import ldata
from Mutual_information import *
from matriz_confusion_presicion import *
from funciones_tp3 import *
#%%
# =============================================================================
# MUTUAL INFORMATION
# =============================================================================
dolph = read_gml('Tp3/dolphins.gml')    
lista = ["infomap","label_prop", "fastgreedy", "eigenvector", "louvain"
        , "edge_betweenness", "walktrap"]

metodos = ["fastgreedy", "eigenvector","edge_betweenness", "louvain", "walktrap", "infomap", "label_prop"]

tabla = np.zeros([len(metodos), len(metodos)])
for i in range(len(metodos)):
    for j in range(len(metodos)):
        tabla[i, j] = round(I_M(dolph, metodos[i], metodos[j]), 4)
import copy
tabla2 = copy.deepcopy(tabla)

tabla2 = (tabla + np.transpose(tabla)) / 2

cuadro = pd.DataFrame(tabla2, columns = metodos, index = metodos)
genders = dict(ldata('Tp3/dolphinsGender.txt'))

# Agrego los sexos a los dicts de cada delfín
for nodo, dict_nodo in dict(dolph.nodes).items():
    dict_nodo['gender'] = genders[nodo] # agrego el sexo del delfín a su dict
#    print('Key = {}, Value = {}'.format(nodo, dict_nodo)) # para chequear que anda

tabla = np.zeros(len(metodos))
for i in range(len(metodos)):
    tabla[i] = I_M(dolph, 'asd', metodos[i], genero = True)
np.save('tabla_info_mutua_generos', tabla)
cuadro = pd.DataFrame(tabla, columns = ['Genero'], index = metodos)
#%%
# =============================================================================
# MATRIZ CONFUSION / MATRIX DE PRESICIONES
# =============================================================================
dolph = read_gml('Tp3/dolphins.gml')    
lista_de_metodos = ["infomap","label_prop", "fastgreedy", "eigenvector",
                    "louvain", "edge_betweenness", "walktrap"] 
part_1 = calcular_particion(dolph, method="fastgreedy", only_labels = True)
part_2 = calcular_particion(dolph, method="walktrap", only_labels = True)
matriz_confusion, presicion = matriz_de_confusion(part_1,part_2,norm=True)
lista_de_metodos_chetos = ["Infomap","Label Propagation", "Fastgreedy",
                           "Eigenvector", "Louvain", "Edge Betweenness",
                           "Walktrap"] 
mp = matriz_de_presiciones(dolph, lista_de_metodos)
plt.figure()
plot_matrix(matriz_confusion, annot=True)
plt.figure()
plot_matrix(mp, lista_de_metodos_chetos, annot=True)
