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
def agregar_centralidad(G, criterio='degree', max_iter=100, tol=1e-3):
    if criterio == 'degree':
        nodes = list(G.nodes())
        degree = list(dict(nx.degree(G)).values())
        return nodes,degree
    if criterio == 'eigen':
        nodes = list(G.nodes())
        eigen = list(dict(nx.eigenvector_centrality(G,
                                                    max_iter=max_iter,
                                                    tol=tol)).values())
        return nodes,eigen    
    if criterio == 'sub':
        nodes = list(G.nodes())
        sub = list(dict(nx.subgraph_centrality(G)).values())
        return nodes,sub 
    if criterio == 'bet':
        nodes = list(G.nodes())
        bet = list(dict(nx.betweenness_centrality(G)).values())
        return nodes,bet
    if criterio == 'flow':
        nodes = list(G.nodes())
        flow = list(dict(nx.current_flow_betweenness_centrality(G)).values())
        return nodes,flow
    if criterio == 'random':
        nodes = list(G.nodes())
        value=random.choice(nodes)
        return value, nodes

def desarme_por_centralidad(g, criterio, parada=2, max_iter=100, tol=1e-3):
    '''
    Toma una red, y bajo algun criterio, elimina nodos de la red iterativamente.
    
    Devuelve 
    --------
    La lista de cantidad de nodos sacados y la lista del tamaño de la componente
    gigante en cada iteracion, ambas normalizadas a la red original.
    
    Sirve para: Degree, Eigenvalues, Subgraph, Betweenness
        
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
    i = 0
    while maxcomp > parada:
        i += 1
#        if i % 10 == 0:
#            print(i, ' pasos')
        nodes, centralidades = agregar_centralidad(G, criterio, max_iter=max_iter)
        j = centralidades.index(max(centralidades))
        G.remove_node(nodes[j])
        centralidades.pop(j)
        nodes.pop(j) #Recorro la lista de nodos eligiendo los que tienen la centralidad máxima            
        #calculo la componente gigante y su tamaño
        cg = max(nx.connected_component_subgraphs(G), key=len)
        maxcomp = cg.order()
        if G.order() > 2:
            cant_nodos_red.append(cant_nodos_red[-1]+(1/nodos_totales))
            cant_nodos_cg.append(maxcomp/cg_original)
        else:
            break
    print(i)
    return cant_nodos_red, cant_nodos_cg
    
def desarme_por_centralidad_flow(g, limite = 1, criterio='flow',parada=2):
    '''
    Toma una red, y bajo algun criterio, elimina nodos de la red iterativamente.
    
    Devuelve 
    --------
    La lista de cantidad de nodos sacados y la lista del tamaño de la componente
    gigante en cada iteracion, ambas normalizadas a la red original.
    
    Sirve para: Flow
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
    i = 0
    while maxcomp > parada:
        i += 1
        G= max(nx.connected_component_subgraphs(G), key=len)
        nodes, degree = agregar_centralidad(G, criterio)         
        for _ in range(limite): #voy a sacar un nodo por iterada
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
    print(i)
    return cant_nodos_red, cant_nodos_cg

def desarme_por_centralidad_aproximado(g, limite = 1, criterio='degree', parada=2):
    '''
    Toma una red, y bajo algun criterio, elimina nodos de la red iterativamente.
    
    Devuelve 
    --------
    La lista de cantidad de nodos sacados y la lista del tamaño de la componente
    gigante en cada iteracion, ambas normalizadas a la red original.
    
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
    i = 0
    while maxcomp > parada:
        i += 1
        nodes, centralidades = agregar_centralidad(G, criterio)
        centralidades, nodes = (list(t) for t in zip(*sorted(zip(centralidades, nodes))))
        for _ in range(limite): #voy a sacar un nodo por iterada
            G.remove_node(nodes[-1])
            centralidades = centralidades[:-1]
            nodes = nodes[:-1]
            #calculo la componente gigante y su tamaño
            maxcomp = max(nx.connected_component_subgraphs(G), key=len).order()
            if G.order() >= 1:
                cant_nodos_red.append(cant_nodos_red[-1]+(1/nodos_totales))
                cant_nodos_cg.append(maxcomp/cg_original)
            else:
                break
    print(i)
    return cant_nodos_red, cant_nodos_cg

#Ojo:Aun no esta probada!!
def desarme_por_centralidad_random(g, limite = 1, criterio='random',parada=2,
                                   n_historias = 100):
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
    nodos_final = np.tile(np.arange(nodos_totales) / nodos_totales, (n_historias, 1))
    cg_final = np.zeros((n_historias, nodos_totales))
    for k in range(n_historias):
        #Copio la red para hacerla bolsa
        G = g.copy()
        G = max(nx.connected_component_subgraphs(G), key=len)
        #Defino las listas que voy a plotear
        cant_nodos_cg = [] #cant de nodos de la cg/cant de nodos de la cg original
        cant_nodos_cg.append(1)
        maxcomp = cg_original
        i = 0
        while maxcomp > parada:
            i += 1
            G= max(nx.connected_component_subgraphs(G), key=len)
            nodo_elegido, nodes = agregar_centralidad(G, criterio)         
            for _ in range(limite): #voy a sacar un nodo por iterada
                j = nodes.index(nodo_elegido)
                G.remove_node(nodes[j])
                nodes.pop(j) #Recorro la lista de nodos eligiendo los que tienen la centralidad máxima            
                #calculo la componente gigante y su tamaño
                cg = max(nx.connected_component_subgraphs(G), key=len)
                maxcomp = cg.order()
                if maxcomp > parada:
                    cant_nodos_cg.append(maxcomp/cg_original)
                else:
                    break
        cg_final[k, :len(cant_nodos_cg)] = cant_nodos_cg
    print(i)
#    import pdb; pdb.set_trace()
    return nodos_final, cg_final