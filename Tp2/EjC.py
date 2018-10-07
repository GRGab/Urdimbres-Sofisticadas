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

#%%
apms = ldata('Tp2/tc02Data/yeast_AP-MS.txt')

y2h = ldata('Tp2/tc02Data/yeast_Y2H.txt')

lit = ldata('Tp2/tc02Data/yeast_LIT.txt')

#lit_r = ldata('Tp2/tc02Data/yeast_LIT_Reguly.txt')


g_apms = nx.Graph()
g_apms.add_edges_from(apms)
#
g_lit = nx.Graph()
g_lit.add_edges_from(lit)

#g_lit_reg = nx.Graph()
#g_lit_reg.add_edges_from(lit_r)

g_y2h = nx.Graph()
g_y2h.add_edges_from(y2h)

#ess = ldata('Tp2/tc02Data/Essential_ORFs_paperHe.txt')

#%%
#    Ejercico c
    
#Criterios
def ordenar_por_centralidad(G, criterio='degree'):
    if criterio == 'degree':
        degree = []
        nodes = []
        for a in G.nodes.keys():
            nodes.append(a)
            degree.append(nx.degree(G)[a])
        return nodes,degree    
    if criterio == 'eigen':
        eigen = []
        nodes = []
        for a in G.nodes.keys():
            nodes.append(a)
            eigen.append(nx.eigenvector_centrality(G)[a])
        return nodes,eigen    
    if criterio == 'sub':
        sub = []
        nodes = []
        for a in G.nodes.keys():
            nodes.append(a)
            sub.append(nx.subgraph_centrality(G)[a])
        return nodes,sub        
    if criterio == 'bet':
        bet = []
        nodes = []
        for a in G.nodes.keys():
            nodes.append(a)
            bet.append(nx.betweenness_centrality(g_apms)[a])
        return nodes, bet   
    if criterio == 'random':
        pass


def desarme_por_grados_1_iteracion(G,limite,criterio='degree',cdp = 1):
    '''
    Toma una red y lista los nodos segun algun criterio de centralidad 
    de meyor a menor (por ejemplo grado de cada nodo). En cada iteracion,
    saco un nodo segun ese criterio. Esto se repite hasta la iteracion 
    'limite'. La funcion devuelve el ultimo grafo, un listado con la cantidad
    acumulativa de nodos sacados de la red y la cantidad de nodos que quedan en la 
    componente gigante.
    '''
    tamaño_cg = [] #lo que va en el eje y
    tamaño_red = [] #lo que va en el eje x
    
    nodes, degree = ordenar_por_centralidad(G, criterio)
    for i in range(1,limite+1): #voy a sacar un nodo por iterada
        j = degree.index(max(degree))
        G.remove_node(nodes[j])
        degree.pop(j)
        nodes.pop(j)
        #calculo la componente gigante y su tamaño
        cg = max(nx.connected_component_subgraphs(G), key=len)
        tamaño_cg.append(cg.order())
        #calculo la cantidad de nodos que saco
        tamaño_red.append(i+cdp*limite)
        if (nx.connected_component_subgraphs(G) == []) :
            print('Me quede sin red.')
            break
        else:
            if (i==limite):
                plt.plot(tamaño_red, tamaño_cg)
                return G, tamaño_red, tamaño_cg 

def desarme_por_grados(g,cant_de_pasos, limite = 5, criterio='degree'):
    '''
    Toma una red, y utlizando la funcion desarme_por_grados_1_iteracion y bajo
    algun criterio, elimina nodos de la red. La variable cant_de_pasos cuenta
    las veces que voy a iterar redefiniendo las centralidades, mientras limite
    refiere a la cantidad de iteradas sin redefinir las centralidades. El total
    de iteraciones es cant_de_pasos*limite. Devuelve la lista de cantidad de
    nodos sacados y la lista del tamaño de la componente gigante en cada 
    iteracion,ambas normalizadas a la red original.
    '''
    #Parametros para la normalizacion
    nodos_totales = g.order()
    cg_original = max(nx.connected_component_subgraphs(g), key=len).order()
    #Copio la red para hacerla bolsa
    G = g.copy()
    #Defino las listas que voy a plotear
    cant_nodos_red = [] #cant de nodos que saque/cant de nodos totales
    cant_nodos_cg = [] #cant de nodos de la cg/cant de nodos de la cg original
    for k in range (0,cant_de_pasos):
        #calculo la fraccion de nodos que quedan y la appendeo
        #cant de nodos que saque
        G_i, nodos_red_i, nodos_cg_i = desarme_por_grados_1_iteracion(
                G, limite, criterio, k)
        for j in range(len(nodos_red_i)):    
            cant_nodos_red.append(nodos_red_i[j]/nodos_totales)
            cant_nodos_cg.append(nodos_cg_i[j]/cg_original)
        G = G_i    
    return cant_nodos_red, cant_nodos_cg
#%%
#    cant_nodos_red, cant_nodos_cg = desarme_por_grados(g_apms,10, criterio = 'degree')
#    plt.plot(cant_nodos_red, cant_nodos_cg)
    cant_nodos_red, cant_nodos_cg = desarme_por_grados(g_apms,15, criterio = 'degree')
    plt.plot(cant_nodos_red, cant_nodos_cg)
    
#%% No darle bola a esto
#ess = ldata('Tp2/tc02Data/Essential_ORFs_paperHe.txt')
#
#ess_1 = []
#
#
#for i in range(2,1158):
#    ess_1.append(ess[i][1]) 
#def asociar_esenciabilidad(red):
#    #Probar esto
#    #lista_de_esenciales = [ list(red.nodes.keys()) for i, j 
#    #                       in list(red.nodes.keys())[j] if 
#    #                       ess_1[i]==list(red.nodes.keys())[j]]
#    lista_de_esenciales = []
#    for i in range(0,len(ess_1)):    
#        for j in range(0,len(list(red.nodes.keys()))):
#            if (ess_1[i]==list(red.nodes.keys())[j]):
#              lista_de_esenciales.append(list(red.nodes.keys())[j])
#    return lista_de_esenciales
#def agregar_importancia(G):
#    g = list(G.nodes())
#    value = np.zeros([len(g)])
#    for h in range(len(g)):
#        for i in range(len(ess)):
#            for j in range(len(ess[i])):
#                if g[h] == ess[i][j]:
#                    value[h] = 1
#                    break
#    return g, value, len(value)
    
    
    
    
    
