#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 10 18:00:08 2018

@author: matias
"""

import networkx as nx
import numpy as np
import random

#Criterios
def ordenar_por_centralidad(G, criterio='degree'):
    if criterio == 'degree':
        nodes = list(G.nodes())
        degree = list(dict(nx.degree(G)).values())  #Calculo la centralidad de mis nodos
        return nodes,degree    
    if criterio == 'eigen':
        nodes = list(G.nodes())
        eigen = list(dict(nx.eigenvector_centrality(G)).values())  #Calculo la centralidad de mis nodos
        return nodes,eigen    
    if criterio == 'sub':
        nodes = list(G.nodes())
        sub = list(dict(nx.subgraph_centrality(G)).values())  #Calculo la centralidad de mis nodos
        return nodes,sub 
    if criterio == 'bet':
        nodes = list(G.nodes())
        bet = list(dict(nx.betweenness_centrality(G)).values())  #Calculo la centralidad de mis nodos
        return nodes,bet
    if criterio == 'flow':
        nodes = list(G.nodes())
        flow = list(dict(nx.current_flow_betweenness_centrality(G)).values())  #Calculo la centralidad de mis nodos
        return nodes,flow
    if criterio == 'random':
        nodes = list(G.nodes())
        value=random.choice(nodes)
        return value, nodes


def desarme_por_grados(g, limite = 1, criterio='degree',parada=2):
    '''
    Toma una red, y utlizando la funcion desarme_por_grados_1_iteracion y bajo
    algun criterio, elimina nodos de la red. La variable cant_de_pasos cuenta
    las veces que voy a iterar redefiniendo las centralidades, mientras limite
    refiere a la cantidad de iteradas sin redefinir las centralidades. El total
    de iteraciones es cant_de_pasos*limite. Devuelve la lista de cantidad de
    nodos sacados y la lista del tamaño de la componente gigante en cada 
    iteracion,ambas normalizadas a la red original.
    Sirve para:Degree, Eigenvalues, Subgraph, Betweenness
        
        '''
    #Parametros para la normalizacion
    nodos_totales = g.order()
    cg_original = max(nx.connected_component_subgraphs(g), key=len).order()
    #Copio la red para hacerla bolsa
    G = g.copy()
    #Defino las listas que voy a plotear
    cant_nodos_red = [] #cant de nodos que saque/cant de nodos totales
    cant_nodos_cg = [] #cant de nodos de la cg/cant de nodos de la cg original
    cant_nodos_red.append(0)
    cant_nodos_cg.append(1)
    maxcomp = cg_original
    while maxcomp > 2:
        nodes, degree = ordenar_por_centralidad(G, criterio)         
        for i in range(1,limite+1): #voy a sacar un nodo por iterada
            j = degree.index(max(degree))
            G.remove_node(nodes[j])
            degree.pop(j)
            nodes.pop(j) #Recorro la lista de nodos eligiendo los que tienen la centralidad máxima            
            #calculo la componente gigante y su tamaño
            cg = max(nx.connected_component_subgraphs(G), key=len)
            maxcomp = cg.order()
            if G.order() > 2:
                cant_nodos_red.append(cant_nodos_red[-1]+(1/nodos_totales))
                cant_nodos_cg.append(maxcomp/cg_original)
            else:
                break
    return cant_nodos_red, cant_nodos_cg
    
def desarme_por_grados_flow(g, limite = 1, criterio='flow',parada=2):
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
    G= max(nx.connected_component_subgraphs(G), key=len)
    #Defino las listas que voy a plotear
    cant_nodos_red = [] #cant de nodos que saque/cant de nodos totales
    cant_nodos_cg = [] #cant de nodos de la cg/cant de nodos de la cg original
    cant_nodos_red.append(0)
    cant_nodos_cg.append(1)
    maxcomp = cg_original
    while maxcomp > parada:
        G= max(nx.connected_component_subgraphs(G), key=len)
        nodes, degree = ordenar_por_centralidad(G, criterio)         
        for i in range(1,limite+1): #voy a sacar un nodo por iterada
            j = degree.index(max(degree))
            G.remove_node(nodes[j])
            degree.pop(j)
            nodes.pop(j) #Recorro la lista de nodos eligiendo los que tienen la centralidad máxima            
            #calculo la componente gigante y su tamaño
            cg = max(nx.connected_component_subgraphs(G), key=len)
            maxcomp = cg.order()
            if maxcomp > parada:
                cant_nodos_red.append(cant_nodos_red[-1]+(1/nodos_totales))
                cant_nodos_cg.append(maxcomp/cg_original)
            else:
                break
    return cant_nodos_red, cant_nodos_cg


#Ojo:Aun no esta probada!!
def desarme_por_grados_random(g, limite = 1, criterio='random',parada=2):
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
    nodos_final = []
    cg_final = []
    for k in range(10):
        #Copio la red para hacerla bolsa
        G = g.copy()
        G= max(nx.connected_component_subgraphs(G), key=len)
        #Defino las listas que voy a plotear
        cant_nodos_red = [] #cant de nodos que saque/cant de nodos totales
        cant_nodos_cg = [] #cant de nodos de la cg/cant de nodos de la cg original
        cant_nodos_red.append(0)
        cant_nodos_cg.append(1)
        maxcomp = cg_original
        while maxcomp > parada:
            G= max(nx.connected_component_subgraphs(G), key=len)
            nodo_elegido, nodes = ordenar_por_centralidad(G, criterio)         
            for i in range(1,limite+1): #voy a sacar un nodo por iterada
                j = nodes.index(nodo_elegido)
                G.remove_node(nodes[j])
                nodes.pop(j) #Recorro la lista de nodos eligiendo los que tienen la centralidad máxima            
                #calculo la componente gigante y su tamaño
                cg = max(nx.connected_component_subgraphs(G), key=len)
                maxcomp = cg.order()
                if maxcomp > parada:
                    cant_nodos_red.append(cant_nodos_red[-1]+(1/nodos_totales))
                    cant_nodos_cg.append(maxcomp/cg_original)
                else:
                    nodos_final.append(cant_nodos_red)
                    cg_final.append(cant_nodos_cg)
                    break
    np.array(nodos_final)
    nodos_final=np.sum(nodos_final,axis=1)/10
    np.array(cg_final)        
    cg_final=np.sum(cg_final,axis=1)/10
    return nodos_final, cg_final