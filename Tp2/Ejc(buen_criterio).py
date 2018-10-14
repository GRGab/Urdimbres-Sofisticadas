#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 13 20:38:20 2018

@author: matias
"""

from lectura import ldata
import networkx as nx
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import random
import collections

import sys
sys.path.append('./Tp2/')
from agregar_esencialidad import agregar_esencialidad
from histograma import histograma
from funciones_de_desarme_ejc import (ordenar_por_centralidad, 
                                        desarme_por_grados,
                                        desarme_por_grados_flow,
                                        desarme_por_grados_random)

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
    
def elimino_nodos(dict_total,k,valores_a_borrar,origen=[]):
    '''Toma el diccionario con los grados como keys y los nodos que le correspoonden como values.
    Si se mete en esta funcion de una (len(origen)=0), quiere decir que la bolsa de nodos
    inicial era suficientemente grande, por lo que simplemente va a eliminar esos nodos de la entrada
    k del diccionario, si eventualmente esa entrada queda vacia, elimina esa entrada del dict.
    Si len(origen)!=0, significa que estan metidos en la bolsa nodos con distinto grado que k.
    Estos nodos tambien hay que eliminarlos del diccionario. La lista con los grados de los nodos
    de la bolsa es origen.'''
    if len(origen)==0:
        dd = list(dict_total[k])
        for u in eliminados:
            dd.remove(u)
        if len(dd)==0:
            dict_total.pop(k)
        else:    
            dict_total[k]=dd                        
    else:                   
        origen.append(k) 
        origen = list(set(origen))
        origen.sort() #ordeno los grados 
        for p in origen:
            dd = list(dict_total[p])
            for l in valores_a_borrar:
                if l in dd:
                    dd.remove(l)
            if len(dd)==0:
                dict_total.pop(p)
                break
            else:
                dict_total[p]=dd
        
    return dict_total


#%%
parada=20
g=g_y2h
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
#import pdb; pdb.set_trace()
listita = list(cant_por_grados_por_sacar.keys())
listita=sorted(listita, reverse=True)
##############
#%%

for kk in listita[:-1]: #sobre todo grado que tengo que sacar
    los_que_saco=cant_por_grados_por_sacar[kk] #es un numerito
    if kk in dict_total.keys(): 
        k=kk
        nodos_para_elegir = dict_total[k] #lista de nodos
#        import pdb; pdb.set_trace()
        i=1 #OJO esto tiene que ser 1, no como en el de abajo
        contador = 0
        conta = 0
        origen= []
        while i!=-1:
#            import pdb; pdb.set_trace()
            if len(nodos_para_elegir) > los_que_saco: #si la bolsa tiene mas elementos de los que necesito
                values = np.random.choice(nodos_para_elegir, size = los_que_saco, replace=False)
                contador+=1
                for n in values:
                    valores_para_sacar.append(n)
                eliminados = list(values)  #Eliminamos del dict los nodos que ya usamos 
                dict_total = elimino_nodos(dict_total, k, origen, eliminados)
                values = []
                i=-1
            else:
                aa = list(dict_total.keys()) #Grados disponibles
                jj = ordenar_grados(aa,k,cercania=i)   
                origen.append(jj)
#                import pdb; pdb.set_trace()
                for m in list(dict_total[jj]):
#                    import pdb; pdb.set_trace()
                    nodos_para_elegir.append(m)
#                    import pdb; pdb.set_trace()        
                i+=1
    else:
        a = list(dict_total.keys()) #Grados disponibles
        j = ordenar_grados(a,kk,cercania=i)       
        nodos_para_elegir = dict_total[j]
        k = j 
        nodos_para_elegir = dict_total[k] #lista de nodos
        i=0 #OJO esto es cero
        conta_1 = 0 
        contador = 0
        origen= []
        while i!=-1:
            if len(nodos_para_elegir) > los_que_saco: #si la bolsa tiene mas elementos de los que necesito
                values = np.random.choice(nodos_para_elegir, size = los_que_saco, replace=False)
                contador+=1
                for n in values:
                    valores_para_sacar.append(n)
                eliminados = list(values)  #Eliminamos del dict los nodos que ya usamos 
                dict_total = elimino_nodos(dict_total, k, origen, eliminados)
                values = []
                i=-1
            else:
                aa = list(dict_total.keys()) #Grados disponibles
                jj = ordenar_grados(aa,k,cercania=i)   
                origen.append(jj)
#                import pdb; pdb.set_trace()
                for m in list(dict_total[jj]):
#                    import pdb; pdb.set_trace()
                    nodos_para_elegir.append(m)
#                    import pdb; pdb.set_trace()                
                i+=1
        
G = g.copy()    
print(conta)
G.remove_nodes_from(valores_para_sacar)
#    orden = G.order()
cg_1 = max(nx.connected_component_subgraphs(G), key=len)
maxcomp_1 = cg_1.order()

b = maxcomp_1/cg_original #Lo que da sacando nodos aleatorios
print(b)
