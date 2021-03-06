#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 28 18:42:00 2018

@author: matias
"""
#Ejercicio 1-c (Presicion/matriz de confusion)

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import networkx as nx
import itertools
import igraph as igraph
from networkx import NetworkXError
from lectura import ldata
from networkx.algorithms.community.community_utils import is_partition
from itertools import product
from networkx.readwrite.gml import read_gml
from histograma import histograma
import seaborn as sns

import sys
sys.path.append('./Tp3/')
from funciones_tp3 import calcular_particion
#%%
def pares(stuff):   
    '''Dada una lista de numeros, devuelve una lista con todos los pares posibles
    a partir de los elementos de la lista. El orden de los pares siempre es 
    el mismo para diferentes listas, esto implica que al comparar distintas
    listas se puede comparar el par de nodos i-esimo en ambas listas
    sin ambiguedad.'''
    lista_de_pares = []
    for i in itertools.combinations(stuff,2):
        lista_de_pares.append(i)
    return lista_de_pares    

def particionar_por_genero_en_lista(G, orden={'f':0, 'NA':1, 'm':2}):
    lista_particiones = []
    for key in dict(G.nodes).keys():
        gender = G.nodes[key]['gender']
        lista_particiones.append(orden[gender])
    return lista_particiones



def matriz_de_confusion(particion_1, particion_2, norm = False):
    '''Dadas dos particiones en forma de listas de numeros, devuelve
    la matriz de confusion correspondiente, segun la definicion vista en clase.
    Se definen las observaciones como los pares de nodos del grafo y las clases
    como los nodos que son del mismo o de diferente tipo.
    (Chequear lo que sique)
    El elemento [i,i] (en la diagonal)
    representa cuantos elementos de la particion de referencia del cluster i se
    encuentran en el cluster i de la particion secundaria. Para el elemento
    de la antidiagonal [i,j], computamos cuantos elementos de la particion de 
    referencia del cluster i no se encuentran en el cluster j. '''

    matriz = np.zeros((2,2))
    pares_a = pares(particion_1)
    pares_b = pares(particion_2)    
    contador_a = 0
    contador_b = 0
    contador_c = 0
    contador_d = 0
    for k in range (len(pares_a)):
        if pares_a[k][0]==pares_a[k][1]:
            if pares_b[k][0]==pares_b[k][1]:
                contador_a +=1
            else:
                contador_c +=1
        if pares_a[k][0]!=pares_a[k][1]:
            if pares_b[k][0]==pares_b[k][1]:
                contador_b +=1
            else:
                contador_d +=1
                
    matriz[0,0] = contador_a
    matriz[0,1] = contador_b
    matriz[1,0] = contador_c
    matriz[1,1] = contador_d 
    presicion = matriz.trace()/matriz.sum()
    
    if norm==True:
        return matriz/matriz.sum(), presicion
    else:
        return matriz, presicion   

def matriz_de_presiciones (graph, lista_de_metodos):    
    '''A partir de un grafo y una lista de metodos a analizar,
    devuelve una matriz en donde el elemento i-j es una comparacion
    entre el metodo i con el j.'''
    diccionario = {} #Claves-->nombres de los metodos. Values-->indice del metodo
    for k, metodo in enumerate(lista_de_metodos):
        diccionario[metodo] = k
    pares_de_metodos = pares(lista_de_metodos)
    matriz_presicion = np.ones((len(lista_de_metodos), len(lista_de_metodos)))
    for i in range(len(pares_de_metodos)):
        part_1 = calcular_particion(graph, method=pares_de_metodos[i][0],
                                    only_labels = True)
        part_2 = calcular_particion(graph, method=pares_de_metodos[i][1],
                                    only_labels = True)
        _, presicion = matriz_de_confusion(part_1,part_2)
        ind_1 = diccionario[pares_de_metodos[i][0]]
        ind_2 = diccionario[pares_de_metodos[i][1]]   
        matriz_presicion[ind_1,ind_2] = presicion
        matriz_presicion[ind_2,ind_1] = presicion    
    return matriz_presicion

def plot_matrix(matriz,
                xlabel=None, ylabel=None, titulo=None,
                labelsize_x=14, labelsize_y=18,
                etiquetas=None, annot = False):
    
    df = pd.DataFrame(matriz, columns = etiquetas, index = etiquetas)
    ax = sns.heatmap(df, annot=annot, fmt=".2f", vmin = 0, vmax=1
                     , square=True)
    ax.set_yticklabels(ax.get_yticklabels(), rotation = 0)
    ax.set_xticklabels(ax.get_xticklabels(), rotation = 0)
#    for label in ax.get_yticklabels():
#            label.set_size(7)
#    for label in ax.get_xticklabels():
#            label.set_size(7)
    for label in ax.get_yticklabels():
        label.set_size(labelsize_y)
    for label in ax.get_xticklabels():
        label.set_size(labelsize_x)
    if titulo is not None:
        ax.set_title(titulo, fontsize=max(labelsize_x, labelsize_y))
    if xlabel is not None:
        ax.set_xlabel(xlabel, fontsize=labelsize_x)
    if ylabel is not None:
        ax.set_ylabel(ylabel, fontsize=labelsize_y)
    plt.tight_layout()
    plt.show()
#%% Pruebita
if __name__ == '__main__':
    a=[1,1,1,2,2,2]
    b=[1,2,1,3,3,3]
    matriz, presicion = matriz_de_confusion(a,b,norm=True)
    plot_matrix(matriz)
#    df_confusion = matriz_de_confusion_wiki(a,b, norm=False)
#    df_confusion = matriz_de_confusion_wiki(part_1,part_2,norm=False)
#    plot_confusion_matrix(df_confusion)
#    print (df_confusion)