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
import collections

from os.path import join as osjoin
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

nombres = {g_apms: "apms", g_lit: "lit", g_lit_reg: "lit_reg",
           g_y2h: "y2h"}
curvas = {}
for g in [g_apms, g_lit, g_lit_reg, g_y2h]:
    for criterio in ['degree', 'eigen', 'sub', 'bet', 'flow', 'random']:
        filename = 'curva_desarme_{}_{}.npz'.format(nombres[g], criterio)
        data = np.load(osjoin('Tp2', 'tc02Data', filename))
        curvas['{}_{}'.format(nombres[g], criterio)] = [data['cant_nodos_red'],
                                                        data['cant_nodos_cg']]
        
#%% Graficar para cada red por separado
        
fontsize = 18
ticksize = 16
        
for g in [g_apms, g_lit, g_lit_reg, g_y2h]:
    with plt.style.context(('seaborn')):
        fig, ax = plt.subplots(figsize=(12,8))
    for criterio in ['degree', 'eigen', 'sub', 'bet', 'flow', 'random']:
        xs = curvas['{}_{}'.format(nombres[g], criterio)][0]
        ys = curvas['{}_{}'.format(nombres[g], criterio)][1]
        if criterio is 'random': # Promedio las 100 historias
            xs = np.average(xs, axis=0)
            ys = np.average(ys, axis=0)
        ax.plot(xs, ys, '.-', label=criterio)
    # Acá falta plotear el puntito de remoción de esenciales
    ax.tick_params(labelsize=ticksize)
    ax.set_xlabel('Fracción de nodos remanentes',
              fontsize=fontsize)
    ax.set_ylabel('Tamaño relativo de la componente gigante',
              fontsize=fontsize)
    ax.legend(fontsize=fontsize)
#    ax.set_xscale('log')
    fig.tight_layout()
#%%
#Punto II
#1015 vs 1179
def desarme_esenciales_saco_ess(g,ess):
    '''
    Toma una red y una lista de nodos esenciales. Primero, elimina los nodos
    esenciales. Luego, elimina nodos al azar pero con el mismo grado
    que los nodos esenciales. En ambos casos se calcula la fraccion de nodos 
    que sobreviven de la componente gigante.
    '''
    ##Primero calculamos la fraccion sacando proteinas esenciales
   
    #Parametros para la normalizacion
    cg_original = max(nx.connected_component_subgraphs(g), key=len).order()
    #Copio la red para hacerla bolsa
    G = g.copy()
    G.remove_nodes_from(ess)
    #calculo la componente gigante y su tamaño
    cg = max(nx.connected_component_subgraphs(G), key=len)
    maxcomp = cg.order()
    a = maxcomp/cg_original #Lo que da sacando esenciales
    return a
#%%
def ordenar_grados(lista,k,cercania=0):
    '''dada una lista de grados, devuelve el grado mas cercano al grado k. 
    Si cercania es n, devuelve el grado en el puesto n de cercania.'''
    from operator import itemgetter

    lista.sort()
    distancia = []
    for i in range (len(lista)):
        distancia.append(abs(lista[i]-k))
    [a,b] = [list(x) for x in zip(*sorted(zip(distancia, lista), key=itemgetter(0)))]    
    return b[cercania]               
                    
                    
                    
                    
def elegir_proximo(lista,m):
    '''Parto de una lista ordenada de menor a mayor, y dado un numero m que no 
    esta en la lista, devuelve el numero de la lista mas cercano a m.'''
    lista_nueva = []
    for p in range (len(lista)):
        lista_nueva.append(abs(lista[p]-m))
    lista_nueva=np.array(lista_nueva)
    ind=list(np.where(lista_nueva==min(lista_nueva))[0])
    ind=ind[0]
    resul = lista[ind]
    return resul               

def desarme_esenciales_saco_random(g,ess):  
    '''
    Luego, elimina nodos al azar pero con el mismo grado
    que los nodos esenciales. En ambos casos se calcula la fraccion de nodos 
    que sobreviven de la componente gigante.
    '''  
    nodos_totales = list(g.nodes()) #lista de nodos de la red
    grados_totales = nx.degree(g) #Dicconario con grados de los nodos de la red
    cg_original = max(nx.connected_component_subgraphs(g), key=len).order()
    
    G_siness = g.copy()
    G_siness.remove_nodes_from(ess)
    grados_totales_siness= nx.degree(G_siness) #Dicconario con grados de los nodos de la red sin ess
        
    lista_grados_esenciales = [] 
    for i in ess:
        if i in nodos_totales:
            a = grados_totales[i]
            lista_grados_esenciales.append(a)
    #Dict con la cant de nodos por sacar con cada grado
    cant_por_grados_por_sacar=dict(collections.Counter(lista_grados_esenciales))  
    
    dict_total = {}
    nodos_no_ess = list(dict(grados_totales_siness).keys())
    grados_no_ess = list(dict(grados_totales_siness).values())
    
    #Me armo un diccionario donde los keys son los grados y los values los nodos que le corresponden
    for i in grados_no_ess:
        nodos_con_k = []
        for j in nodos_no_ess:
            if grados_totales_siness[j] == i:
                nodos_con_k.append(j)
            dict_total[i] = nodos_con_k
    valores_para_sacar = []
    
    ##############
    for kk in cant_por_grados_por_sacar.keys(): #sobre todo grado que tengo que sacar
        los_que_saco=cant_por_grados_por_sacar[kk] #es un numerito
        if kk in dict_total.keys(): 
            k=kk
            nodos_para_elegir = dict_total[k] #lista de nodos
            i=0
            while i!=-1:
                if len(nodos_para_elegir) > los_que_saco: #si la bolsa tiene mas elementos de los que necesito
                    values = np.random.choice(nodos_para_elegir, size = los_que_saco, replace=False)
                    for i in values:
                        valores_para_sacar.append(i)
                    eliminados = list(values)  #Eliminamos del dict los nodos que ya usamos 
                    for l in eliminados:
                        nodos_para_elegir.remove(l)
                    dict_total[k]=nodos_para_elegir
                    i=-1
                else:
                    a = list(dict_total.keys()) #Grados disponibles
                    j = ordenar_grados(a,k,cercania=i)      
                    nodos_para_elegir = dict_total[j]
                    i+=1
                    k=j
        else:
            a = list(dict_total.keys()) #Grados disponibles
            j = ordenar_grados(a,k,cercania=i)       
            nodos_para_elegir = dict_total[j]
            k = j 
            nodos_para_elegir = dict_total[k] #lista de nodos
            i=0
            while i!=-1:
                if len(nodos_para_elegir) > los_que_saco: #si la bolsa tiene mas elementos de los que necesito
                    values = np.random.choice(nodos_para_elegir, size = los_que_saco, replace=False)
                    for i in values:
                        valores_para_sacar.append(i)
                    eliminados = list(values)  #Eliminamos del dict los nodos que ya usamos 
                    for l in eliminados:
                        nodos_para_elegir.remove(l)
                    dict_total[k]=nodos_para_elegir
                    i=-1
                else:
                    a = list(dict_total.keys()) #Grados disponibles
                    j = ordenar_grados(a,k,cercania=i)      
                    nodos_para_elegir = dict_total[j]
                    i+=1
                    k=j
                
    G = g.copy()    
    G.remove_nodes_from(valores_para_sacar)
#    orden = G.order()
    cg_1 = max(nx.connected_component_subgraphs(G), key=len)
    maxcomp_1 = cg_1.order()
    b = maxcomp_1/cg_original #Lo que da sacando nodos aleatorios
    return b

def analisis_desarme_esenciales(g, ess, numero_de_tiradas):
    lista = []
    valor_real = desarme_esenciales_saco_ess(g,ess)
    for i in range(numero_de_tiradas):
        b = desarme_esenciales_saco_random(g,ess)
        lista.append(b)
        
    return lista, valor_real
#%%    
nombre_archivo = 'Tp2/Ej c/lit(1000 tiradas).npz'
import time; ti = time.time()
lista, valor_real = analisis_desarme_esenciales(g_lit, ess, int(1e3))
np.savez(nombre_archivo, lista=lista,
         valor_real=valor_real)
tf = time.time(); print(tf-ti, 'segundos')
#%%
nombre_archivo = 'Tp2/Ej c/lit_reg_(1000 tiradas).npz'
import time; ti = time.time()
lista, valor_real = analisis_desarme_esenciales(g_lit_reg, ess, int(1e3))
np.savez(nombre_archivo, lista=lista,
         valor_real=valor_real)
tf = time.time(); print(tf-ti, 'segundos')
#%%
nombre_archivo = 'Tp2/Ej c/y2h_(1000 tiradas).npz'
import time; ti = time.time()
lista, valor_real = analisis_desarme_esenciales(g_y2h, ess, int(1e3))
np.savez(nombre_archivo, lista=lista,
         valor_real=valor_real)
tf = time.time(); print(tf-ti, 'segundos')
#%%
nombre_archivo = 'Tp2/Ej c/apms_(1000 tiradas).npz'
import time; ti = time.time()
lista, valor_real = analisis_desarme_esenciales(g_apms, ess, int(1e3))
np.savez(nombre_archivo, lista=lista,
         valor_real=valor_real)
tf = time.time(); print(tf-ti, 'segundos')
#%% Cargar datos
nombre_archivo_cargado = 'Tp2/Ej c/apms_nuevo(1000 tiradas).npz'
#nombre_archivo_cargado = 'apms_nuevo(1000 tiradas).npz'
lista = np.load(nombre_archivo_cargado)['lista']
valor_real = np.load(nombre_archivo_cargado)['valor_real']
#%% Ploteamos resultados
fig, ax = histograma(lista, bins=10, density=True, errbars=False, 
                     titulo=r'Distribución de la fraccion de la CG que sobrevive bajo $H_0$ (LIT_REG)',
                     xlabel='Fraccion de la componente gigante')
ax.axvline(valor_real, color='deeppink',
           label='Valor real = {}'.format(valor_real))
ax.legend()
fig.savefig('Tp2/Ej c/lit_reg_1000_tiradas).png', bbox_inches='tight')
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
