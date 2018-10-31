#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 28 18:37:39 2018

@author: matias
"""

import numpy as np
from networkx.readwrite.gml import read_gml
from funciones_tp3 import(guardar_particiones, graficar_dist)
#%%
dolph = read_gml('Tp3/dolphins.gml')    
lista_de_metodos = ["infomap","label_prop", "fastgreedy", "eigenvector",
                    "louvain", "edge_betweenness", "walktrap"]
guardar_particiones(dolph, 1000,200, lista_de_metodos)
#%% Importamos todo y graficamos
npzfile = np.load('Tp3/tc03Data/Ej_b_particiones_numeros.npz')
rewire = npzfile['mod_rewire']
original = npzfile['mod_original']
for i in range(0,7):
    graficar_dist(dolph, lista_de_metodos, rewire[i], original[i],
                  metodo = i)
#%% Pueba de las funciones del punto 1-b
dolph = read_gml('Tp3/dolphins.gml')    
lista_de_metodos = ["infomap","label_prop", "fastgreedy", "eigenvector",
                    "louvain", "edge_betweenness", "walktrap"]
#    guardar_particiones(dolph, 200,10, lista_de_metodos)
