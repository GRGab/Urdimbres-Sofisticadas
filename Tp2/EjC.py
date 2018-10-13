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
import time

import sys
sys.path.append('./Tp2/')
from agregar_esencialidad import agregar_esencialidad
from histograma import histograma
from funciones_de_desarme_ejc import (agregar_centralidad, 
                                        desarme_por_centralidad,
                                        desarme_por_centralidad_flow,
                                        desarme_por_centralidad_random)
#%%
apms = ldata('Tp2/tc02Data/yeast_AP-MS.txt')

y2h = ldata('Tp2/tc02Data/yeast_Y2H.txt')

lit = ldata('Tp2/tc02Data/yeast_LIT.txt')

lit_r = ldata('Tp2/tc02Data/yeast_LIT_Reguly.txt')
lit_r = [fila[:2] for fila in lit_r[1:]]

g_apms = nx.Graph()
g_apms.add_edges_from(apms)

g_lit = nx.Graph()
g_lit.add_edges_from(lit)

g_lit_reg = nx.Graph()
g_lit_reg.add_edges_from(lit_r)

g_y2h = nx.Graph()
g_y2h.add_edges_from(y2h)

ess = ldata('Tp2/tc02Data/Essential_ORFs_paperHe.txt')
ess =  ess[2:-4]
ess = [fila[1] for fila in ess]
ess = np.unique(ess)
#%%
#    Ejercico c
#Punto I
ti = time.time()
cant_nodos_red, cant_nodos_cg = desarme_por_centralidad(g_apms,
                                                        criterio = 'eigen',
                                                        parada=2)
tf = time.time(); print(tf-ti, 'segundos') 
plt.plot(cant_nodos_red, cant_nodos_cg)
#%% Por SUBGRAPH

### APMS
ti = time.time() 
cant_nodos_red, cant_nodos_cg = desarme_por_centralidad(g_apms,
                                                        criterio = 'sub')
tf = time.time(); print(tf-ti, 'segundos') 

plt.plot(cant_nodos_red, cant_nodos_cg)
np.savez('curva_desarme_apms_sub.npz', cant_nodos_red=cant_nodos_red,
         cant_nodos_cg=cant_nodos_cg)

### LIT
ti = time.time() 
cant_nodos_red, cant_nodos_cg = desarme_por_centralidad(g_lit,
                                                        criterio = 'sub')
tf = time.time(); print(tf-ti, 'segundos')
plt.plot(cant_nodos_red, cant_nodos_cg)
np.savez('curva_desarme_lit_sub.npz', cant_nodos_red=cant_nodos_red,
         cant_nodos_cg=cant_nodos_cg)
#%%
### Por grado
ti = time.time()
cant_nodos_red, cant_nodos_cg = desarme_por_centralidad(g_apms,
                                                        criterio = 'degree')
tf = time.time(); print(tf-ti, 'segundos') 
plt.figure(); plt.plot(cant_nodos_red, cant_nodos_cg)
#np.savez('curva_desarme_apms_deg.npz', cant_nodos_red=cant_nodos_red,
#         cant_nodos_cg=cant_nodos_cg)
#%% flow
ti = time.time()
cant_nodos_red, cant_nodos_cg = desarme_por_centralidad_flow(g_apms)
plt.plot(cant_nodos_red, cant_nodos_cg)
tf = time.time(); print(tf-ti, 'segundos')
# 793 segundos
np.savez('curva_desarme_apms_flow.npz', cant_nodos_red=cant_nodos_red,
         cant_nodos_cg=cant_nodos_cg)
#%% betweenness
ti = time.time()
cant_nodos_red, cant_nodos_cg = desarme_por_centralidad_flow(g_apms)
plt.plot(cant_nodos_red, cant_nodos_cg)
tf = time.time(); print(tf-ti, 'segundos')
np.savez('curva_desarme_apms_bet.npz', cant_nodos_red=cant_nodos_red,
         cant_nodos_cg=cant_nodos_cg)

#%% Random

### APMS
ti = time.time()
cant_nodos_red, cant_nodos_cg = desarme_por_centralidad_random(g_apms,
                                                               n_historias=1,
                                                               parada=2)
tf = time.time(); print(tf-ti, 'segundos')
#np.savez('curva_desarme_apms_random.npz', cant_nodos_red=cant_nodos_red,
#         cant_nodos_cg=cant_nodos_cg)
plt.figure()
plt.plot(np.average(cant_nodos_red, axis=0), np.average(cant_nodos_cg, axis=0))
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
    nodos_totales = list(g.nodes()) #lista de nodos de la red
    grados_totales = nx.degree(g) #Dicconario con grados de los nodos de la red
    cg_original = max(nx.connected_component_subgraphs(g), key=len).order()
    
    #Copio la red para hacerla bolsa
    G = g.copy()
    G.remove_nodes_from(ess)
    #calculo la componente gigante y su tamaño
    cg = max(nx.connected_component_subgraphs(G), key=len)
    maxcomp = cg.order()
##########
    
    #Ahora calculamos la fraccion sacando proteinas random con el mismo grado.
    
    GGG = g.copy()
    GGG.remove_nodes_from(ess) ##Extraigo los nodos esenciales
    nodes = list(GGG.nodes()) #lista de nodos de la red no esenciales  
    #Obs, tengo que sacar nodos no esenciales, pero los tengo que sacar de la red original!
    
    GG = g.copy()    
    lista_para_eliminar = []
    for nodito in ess: #itero sobre los nodos esenciales
        if nodito in nodos_totales: #si el nodo esencial esta en la red original
            gradito = grados_totales[nodito]      
            grados = list(dict(nx.degree(GGG)).values())
            grados = np.array(grados)
            j = np.where(grados==gradito)[0] #j son los indices (de GGG) de los nodos con igual grado
                                             #que los nodos esenciales 
            if len(j)>0:
                #poner un print para ver cuando paso eso
                # no puedo agarrar este pruebo el grado de al lado
                lista_de_nodos = []
                for k in j:
                    lista_de_nodos.append(nodes[k])
                    value=random.choice(lista_de_nodos)
                    lista_para_eliminar.append(value)
    GG.remove_nodes_from(lista_para_eliminar)
    cg_1 = max(nx.connected_component_subgraphs(GG), key=len)
    maxcomp_1 = cg_1.order()
    
    a = maxcomp/cg_original #Lo que da sacando esenciales
    b = maxcomp_1/cg_original #Lo que da sacando nodos aleatorios
    return a,b

def analisis_desarme_esenciales(g, ess, numero_de_tiradas):
    lista = []
    valor_real, _ = desarme_esenciales(g,ess)
    for i in range(numero_de_tiradas):
        _, b = desarme_esenciales(g,ess)
        lista.append(b)
    return lista, valor_real

#%%    
nombre_archivo = 'lit_(1000 tiradas).npz'
import time; ti = time.time()
lista, valor_real = analisis_desarme_esenciales(g_lit, ess, int(1e3))
np.savez(nombre_archivo, lista=lista,
         valor_real=valor_real)
tf = time.time(); print(tf-ti, 'segundos')
#%%
nombre_archivo = 'lit_reg_(1000 tiradas).npz'
import time; ti = time.time()
lista, valor_real = analisis_desarme_esenciales(g_lit_reg, ess, int(1e3))
np.savez(nombre_archivo, lista=lista,
         valor_real=valor_real)
tf = time.time(); print(tf-ti, 'segundos')
#%%
nombre_archivo = 'y2h_(1000 tiradas).npz'
import time; ti = time.time()
lista, valor_real = analisis_desarme_esenciales(g_y2h, ess, int(1e3))
np.savez(nombre_archivo, lista=lista,
         valor_real=valor_real)
tf = time.time(); print(tf-ti, 'segundos')
#%%
nombre_archivo = 'apms_(1000 tiradas).npz'
import time; ti = time.time()
lista, valor_real = analisis_desarme_esenciales(g_apms, ess, int(1e3))
np.savez(nombre_archivo, lista=lista,
         valor_real=valor_real)
tf = time.time(); print(tf-ti, 'segundos')
#%% Cargar datos
nombre_archivo_cargado = 'Ej c/lit_reg_(100 tiradas).npz'
lista = np.load(nombre_archivo_cargado)['lista']
valor_real = np.load(nombre_archivo_cargado)['valor_real']
#%% Ploteamos resultados
fig, ax = histograma(lista, bins=5, density=True, errbars=False, 
                     titulo=r'Distribución de la fraccion de la CG que sobrevive bajo $H_0$',
                     xlabel='Fraccion de la componente gigante')
ax.axvline(valor_real, color='deeppink',
           label='Valor real = {}'.format(valor_real))
ax.legend()
#fig.savefig('Ej c/Ej c punto 2 (1000 tiradas).pdf', bbox_inches='tight')
plt.show()
#La esencialidad de los nodos de la red no tienen que ver con el grado que tienen.
#%%

#from datetime import datetime
#datestring = datetime.strftime(datetime.now(), '%Y/%m/%d_%H:%M:%S')

# Guardar archivo
#thefile = open('Ej c/lit_(100 tiradas)_'+datestring+'.txt', 'w+')
#for item in lista:
#  thefile.write("%s\n" % item)
#thefile.close()




#lista = np.loadtxt('Ej c/Ej c punto 2 (1000 tiradas)y2h.txt')
#fig, ax = histograma(lista, bins=10, density=True, errbars=False, 
#                     titulo=r'Distribución de la fraccion de la CG que sobrevive bajo $H_0$',
#                     xlabel='Fraccion de la componente gigante')
