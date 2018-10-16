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

import sys
from time import time
sys.path.append('./Tp2/')
from agregar_esencialidad import agregar_esencialidad
from collections import Counter
from operator import itemgetter

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

def desarme_esenciales_saco_ess(g, ess, devolver_num_nodos_elim=False):
    '''
    Toma una red y una lista de nodos esenciales. Elimina los nodos esenciales
    y calcula el tamaño relativo de la componente gigante en la nueva red.
    '''
   #Parametros para la normalizacion
    cg_original = max(nx.connected_component_subgraphs(g), key=len).order()
    #Copio la red para hacerla bolsa
    G = g.copy()
    G.remove_nodes_from(ess)
    #calculo la componente gigante y su tamaño
    cg = max(nx.connected_component_subgraphs(G), key=len)
    maxcomp = cg.order()
    a = maxcomp/cg_original #Lo que da sacando esenciales
    if devolver_num_nodos_elim:
        num_nodos_elim = g.order() - G.order()
        return a, num_nodos_elim
    else:
        return a

def obtener_mascercano(lista, k, cercania=0):
    '''dada una lista de grados, devuelve el grado mas cercano al grado k. 
    Si cercania es n, devuelve el grado en el puesto n de cercania.'''
    lista.sort()
    distancia = []
    for i in range (len(lista)):
        distancia.append(abs(lista[i]-k))
    [a,b] = [list(x) for x in zip(*sorted(zip(distancia, lista), key=itemgetter(0)))]
    # Siempre devuelve con menor cercanía el grado más próximo por abajo, 
    # y con mayor cercanía el más próximo por arriba
    return b[cercania]


def elegir_proximo(lista,m): # Sería más bien "elegir_anterior"
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
    
def elimino_nodos(dict_total, k, valores_a_borrar, origen=[]):
    '''Toma el diccionario con los grados como keys y los nodos que le correspoonden como values.
    Si se mete en esta funcion de una (len(origen)=0), quiere decir que la bolsa de nodos
    inicial era suficientemente grande, por lo que simplemente va a eliminar esos nodos de la entrada
    k del diccionario, si eventualmente esa entrada queda vacia, elimina esa entrada del dict.
    Si len(origen)!=0, significa que estan metidos en la bolsa nodos con distinto grado que k.
    Estos nodos tambien hay que eliminarlos del diccionario. La lista con los grados de los nodos
    de la bolsa es origen.'''
    
    if len(origen)==0:
        dd = list(dict_total[k])
        for u in valores_a_borrar:
            dd.remove(u)
        if len(dd)==0:
            dict_total.pop(k)
        else:    
            dict_total[k]=dd
    else:                   
        origen.append(k)
        origen = list(set(origen)) # Chequear que no haya repetidos
        origen.sort() #ordeno los grados 
        for p in origen:
            dd = list(dict_total[p])
            for l in valores_a_borrar:
                if l in dd:
                    dd.remove(l)
            if len(dd)==0:
                dict_total.pop(p)
            else:
                dict_total[p]=dd
        
    return dict_total


#%%
def efecto_nodos_equivalentes(g, ess, devolver_num_nodos_elim=False):
    """Elimina un conjunto aleatorio de nodos del grafo 'g' con igual secuencia
    de grados que la del conjunto de nodos esenciales 'ess' y devuelve la fracción
    que representa la componente gigante final respecto de la componente gigante
    inicial.
    
    El resultado es aleatorio, por lo cual es recomendable promediar sobre
    varias historias.
    
    La diferencia con el caso anterior es que ahora vamos a adoptar un criterio
    ligeramente distinto para lidiar con el caso (muy común) en el que no hay
    suficientes nodos del grado deseado para eliminar."""   
        
    grados_esenciales = [] # contiene repetidos
    for nodo in ess:
        if nodo in g.nodes():
            k = g.degree[nodo]
            grados_esenciales.append(k)
    # Dict con la cant de nodos esenciales que tienen cierto grado
    # Esto nos dice cuántos nodos no esenciales hay que eliminar para cada grado
    grado_to_numesenciales = dict(Counter(grados_esenciales))
    # Lista de grados esenciales sin repetir de mayor a menor
    # Sobre esto vamos a iterar para remover los nodos que correspondan
    grados_esenciales = sorted(np.unique(grados_esenciales), reverse=True)
    

    
    # Dict que manda un valor de grado al conjunto de no esenciales con dicho grado
    # De acá vamos sacando nodos y metiéndolos en la bolsa de nodos a eliminar
    nodos_no_esenciales = [nodo for nodo in g.nodes() if nodo not in ess]
#    import pdb; pdb.set_trace()
    # Grados presentes en la red
    grados_presentes = np.unique(list(dict(g.degree()).values()))
    # Inicializamos el diccionario a completar
    grado_to_noesenciales = {}
    # Completamos
    for k in grados_presentes:
        noesenciales_con_grado_k = []
        for nodo in nodos_no_esenciales:
            if g.degree[nodo] == k:
                noesenciales_con_grado_k.append(nodo)
            if len(noesenciales_con_grado_k) != 0:
                grado_to_noesenciales[k] = noesenciales_con_grado_k
#    import pdb; pdb.set_trace()
    # OJO: grado_to_noesenciales es mutable y va cambiando a medida que avanza
    # la ejecución.
    
# =============================================================================
#           CÓDIGO EN CUARENTENA
#
#     G_siness = g.copy()
#     G_siness.remove_nodes_from(ess)
#     dict_grados_siness = dict(nx.degree(G_siness)) #Dicconario con grados de los nodos de la red sin ess
#     grado_to_noesenciales = {}
#     nodos_no_ess = dict_grados_siness.keys()
#     grados_no_ess = dict_grados_siness.values()
#     for k in np.unique(list(grados_no_ess)):
#         nodos_con_grado_k = []
#         for nodo in nodos_no_ess:
#             if dict_grados_siness[nodo] == k:
#                 nodos_con_grado_k.append(nodo)
#             grado_to_noesenciales[k] = nodos_con_grado_k
# =============================================================================
    
    
    nodos_a_eliminar = []
    for k in grados_esenciales: # sobre todo grado del cual tengo que sacar nodos

        # cuantos_sacar es cuántos nodos tengo que sacar en esta iteración
        cuantos_sacar = grado_to_numesenciales[k]
        # Grados disponibles para ser eliminados
        grados_disponibles = list(grado_to_noesenciales.keys())
        # Si no tenemos nodos no esenciales con el grado deseado, saltamos
        # directamente al grado disponible más cercano
        if k not in grados_disponibles: 
            k = obtener_mascercano(grados_disponibles, k, cercania=0) # Grado más cercano

# =============================================================================
#         if k in grados_disponibles: # Si tenemos nodos no esenciales con el grado deseado
#             i = 1 #OJO esto tiene que ser 1, no como en el de abajo
#         else:
#             k = obtener_mascercano(grados_disponibles, k, cercania=0) # Grado más cercano
#             i = 0 #OJO esto es cero
# =============================================================================
        nodos_para_elegir = grado_to_noesenciales[k] # lista de nodos no ess con el grado deseado
        origen = []
        
        # Parámetros del loop while
        procesando_k = True # Cuando se vuelve False, pasamos al siguiente k
        i = 1 # Aumenta cada vez que repito el loop
        while procesando_k:
            # Si la bolsa tiene suficientes elementos
            if len(nodos_para_elegir) >= cuantos_sacar:
                nominados = np.random.choice(nodos_para_elegir, size = cuantos_sacar, replace=False)
                for n in nominados:
                    nodos_a_eliminar.append(n)
                #Eliminamos del dict los nodos que ya usamos:
                eliminados = list(nominados)  
                grado_to_noesenciales = elimino_nodos(grado_to_noesenciales, k, eliminados, origen)
                procesando_k = False # Termina el while, pasamos al siguiente k

            else: 
                # si no hay suficientes nodos no esenciales de grado k, voy a buscar
                # más y los agrego a la bolsa
                k_agregado = obtener_mascercano(grados_disponibles, k, cercania=i)
                origen.append(k_agregado)
                for m in list(grado_to_noesenciales[k_agregado]):
                    nodos_para_elegir.append(m)
                i += 1

    G = g.copy()
    G.remove_nodes_from(nodos_a_eliminar)
    cg_1 = max(nx.connected_component_subgraphs(G), key=len)
    maxcomp_1 = cg_1.order()
    tamaño_cg_original = max(nx.connected_component_subgraphs(g), key=len).order()
    b = maxcomp_1 / tamaño_cg_original #Lo que da sacando nodos aleatorios
    if devolver_num_nodos_elim:
        return b, len(nodos_a_eliminar)
    else:
        return b

#%%
grafos = [g_lit, g_lit_reg, g_apms, g_y2h] # Respetar este orden
nombres = {g_apms: "APMS", g_lit: "Lit", g_lit_reg: "Lit_reg",
           g_y2h: "Y2H"}

# Genero las fracciones correspondientes a remoción de esenciales
fracs_esen = np.zeros((4))
for i, g in enumerate(grafos):
    print('Eliminando nodos esenciales de la red ', nombres[g])
    fracs_esen[i], numnods_elim = desarme_esenciales_saco_ess(g, ess, devolver_num_nodos_elim=True)
    print('# Nodos eliminados: ', numnods_elim)

# Genero las fracciones correspondientes a remoción de equivalentes en grado
fracs_equiv, sigmas = np.zeros((2, 4))
n_historias = 10
ti = time()
for i, g in enumerate(grafos):
    print('Eliminando nodos no esenciales equivalentes de la red ', nombres[g])
    print('Número de historias: ', n_historias)
    resultados = np.zeros((n_historias))
    numnods_elim = np.zeros((n_historias))
    for j in range(n_historias):
        resultados[j], numnods_elim[j] = efecto_nodos_equivalentes(g, ess, devolver_num_nodos_elim=True)
    fracs_equiv[i] = np.average(resultados)
    sigmas[i] = np.std(resultados)
    print('# Nodos eliminados: {} +/- {}'.format(np.average(numnods_elim), np.std(numnods_elim)))
tf = time(); print(tf - ti, ' segundos')

zscores = (fracs_esen - fracs_equiv) / sigmas

#np.savez('Tp2/tc02Data/unosdatos.npz', fracs_esen=fracs_esen,
#         fracs_equiv=fracs_equiv, sigmas=sigmas, zscores=zscores)

# Tabla

equiv_prolijo = ['{:.2f} +/- {:.2f}'.format(x, y) for x, y in zip(fracs_equiv, sigmas)]
pd.set_option('precision', 2)
tabla3 = pd.DataFrame(data={'Esenciales':fracs_esen,
                            'No esenciales equivalentes:':equiv_prolijo,
                            'Z-scores':zscores},
                      index=['Lit', 'Lit_reg', 'APMS', 'Y2H'])
tabla3