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
import random
from histograma import histograma

import sys
sys.path.append('./Tp2/')
from agregar_esencialidad import agregar_esencialidad
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

ess = ldata('Tp2/tc02Data/Essential_ORFs_paperHe.txt')
ess = [fila for fila in ess if (fila != [] and fila[0] != 'updated')]
ess = [fila[1] for fila in ess]
#%%
#    Ejercico c
    
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

#%%

#    cant_nodos_red, cant_nodos_cg = desarme_por_grados(g_apms, limite = 1,criterio = 'eigen', parada=1200)
#    plt.plot(cant_nodos_red, cant_nodos_cg)
#    
#    cant_nodos_red, cant_nodos_cg = desarme_por_grados(g_apms, criterio = 'sub')
#    plt.plot(cant_nodos_red, cant_nodos_cg)
####    
#    cant_nodos_red, cant_nodos_cg = desarme_por_grados(g_apms, criterio = 'degree',parada=1000)
#    plt.plot(cant_nodos_red, cant_nodos_cg)
    
#cant_nodos_red, cant_nodos_cg = desarme_por_grados_flow(g_apms,parada=800)
#plt.plot(cant_nodos_red, cant_nodos_cg)
#    
#%%
#Punto II
def desarme_esenciales(g,ess):
    '''
    Toma una red y una lista de nodos esenciales. Primero, elimina los nodos
    esenciales. Luego, elimina nodos al azar pero con el mismo grado
    que los nodos esenciales. En ambos casos se calcula la fraccion de nodos 
    que sobreviven de la componente gigante.
    '''
    ##Primero calculamos la fraccion sacando proteinas esenciales
   
    #Parametros para la normalizacion
#    nodos_totales = g.order()
    cg_original = max(nx.connected_component_subgraphs(g), key=len).order()
    
    #Copio la red para hacerla bolsa
    G = g.copy()
    G.remove_nodes_from(ess)
    #calculo la componente gigante y su tamaño
    cg = max(nx.connected_component_subgraphs(G), key=len)
    maxcomp = cg.order()
##########
    #Ahora calculamos la fraccion sacando proteinas random con el mismo grado.
    GG = g.copy()
    for nodito in ess: #itero sobre los nodos esenciales
        nodes = list(GG.nodes()) #lista de nodos de la red
        if nodito in list(GG.nodes()):
            grados = nx.degree(GG)
            gradito = grados[nodito]
            
            grados = np.array(grados)
            j = np.where(grados==gradito)[0] #j son los indices de los nodos con igual grado
                                          #que los nodos esenciales 
            lista_de_nodos = []
            for k in j:
                lista_de_nodos.append(nodes[k])
            value=random.choice(lista_de_nodos)
            GG.remove_node(value)
    
    cg_1 = max(nx.connected_component_subgraphs(GG), key=len)
    maxcomp_1 = cg_1.order()
    
    a = maxcomp/cg_original #Lo que da sacando esenciales
    b = maxcomp_1/cg_original #Lo que da sacando nodos aleatorios
    return a,b

def analisis_desarme_esenciales(g, ess):
    lista = []
    valor_real, _ = desarme_esenciales(g,ess)
    for i in range(100):
        _, b = desarme_esenciales(g,ess)
        lista.append(b)
    return lista
#%%
import time; ti = time.time()
lista = analisis_desarme_esenciales(g_apms, ess)
tf = time.time(); print(tf-ti, 'segundos')
#%%
valor_real = 0.3237051792828685
fig, ax = histograma(lista, bins=20, density=True, errbars=False, 
                     titulo=r'Distribución de la fraccion de la cg que sobrevive bajo $H_0$',
                     xlabel='Fraccion de la componente gigante')
ax.axvline(valor_real, color='deeppink',
           label='Valor real = {}'.format(valor_real))
ax.legend()
plt.show()
#La esencialidad de los nodos de la red no tienen que ver con el grado que tienen.


