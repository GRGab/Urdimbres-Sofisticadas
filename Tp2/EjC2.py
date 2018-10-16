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
def efecto_nodos_equivalentes_criterioviejo(g, ess, devolver_num_nodos_elim=False):
    """Elimina un conjunto aleatorio de nodos del grafo 'g' con igual secuencia
    de grados que la del conjunto de nodos esenciales 'ess' y devuelve la fracción
    que representa la componente gigante final respecto de la componente gigante
    inicial.
    
    El resultado es aleatorio, por lo cual es recomendable promediar sobre
    varias historias."""   
    
    # Creamos lista de grados de los nodos esenciales (contiene repetidos)
    grados_esenciales = [] 
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
    # OJO: grado_to_noesenciales es mutable y va cambiando a medida que avanza
    # la ejecución.

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
        # lista de nodos no esenciales con grado k
        nodos_para_elegir = grado_to_noesenciales[k] 
        origen = []
        i = 1 # Aumenta cada vez que repito el loop
        while len(nodos_para_elegir) < cuantos_sacar:
            # si no hay suficientes nodos no esenciales de grado k, voy a buscar
            # más y los agrego a la bolsa
            k_agregado = obtener_mascercano(grados_disponibles, k, cercania=i)
            # Dejamos registrado que k_agregado es uno de los grados de los
            # cuales hay que eliminar los nodos:
            origen.append(k_agregado)
            # Agregamos todos los nodos de grado  k_agregado a la bolsa
            # de nodos para elegir
            nodos_para_elegir += grado_to_noesenciales[k_agregado]
            # Aumentamos i para ir a buscar a grados más lejanos al original,
            # en caso de ser encesario
            i += 1
        
        # Si terminó el while es porque ya hay suficientes no esenciales
        # en la bolsa.
        # Metemos la mano en la bosla y sacamos los nodos a eliminar
        nominados = np.random.choice(nodos_para_elegir, size = cuantos_sacar, replace=False)
        nodos_a_eliminar += list(nominados)
        #Eliminamos del dict los nodos que ya marcamos para eliminar:
        grado_to_noesenciales = elimino_nodos(grado_to_noesenciales, k, nominados, origen)

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
        
def efecto_nodos_equivalentes(g, ess, devolver_num_nodos_elim=False):
    """Elimina un conjunto aleatorio de nodos del grafo 'g' con igual secuencia
    de grados que la del conjunto de nodos esenciales 'ess' y devuelve la fracción
    que representa la componente gigante final respecto de la componente gigante
    inicial.
    
    El resultado es aleatorio, por lo cual es recomendable promediar sobre
    varias historias."""   
    
    # Creamos lista de grados de los nodos esenciales (contiene repetidos)
    grados_esenciales = [] 
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
    # OJO: grado_to_noesenciales es mutable y va cambiando a medida que avanza
    # la ejecución.

    nodos_a_eliminar = []
    for k in grados_esenciales: # sobre todo grado del cual tengo que sacar nodos      
        # Cuántos falta marcar como eliminados
        cuantos_faltan = grado_to_numesenciales[k] 
        # Grados disponibles para eliminar
        grados_disponibles = list(grado_to_noesenciales.keys())     
        if grados_disponibles == []:
            # Si todavía nos faltaban grados esenciales pero ya no queda
            # ningún grado disponible, queremos terminar el for loop
            # instantáneamente
            break
        if k not in grados_disponibles:
            # Si no tenemos nodos no esenciales con el grado deseado, saltamos
            # directamente al grado disponible más cercano
            k = obtener_mascercano(grados_disponibles, k, cercania=0)
        nodos_a_eliminar_temp = []
        while cuantos_faltan > 0:
            n_k = len(grado_to_noesenciales[k])
            if  n_k < cuantos_faltan:
                # Si no hay suficientes no esenciales de grado k,
                # marco todos ellos para ser eliminados (con certeza) y me fijo
                # en los nodos de grado disponible más cercano para completar
                nodos_a_eliminar_temp += grado_to_noesenciales[k]
                cuantos_faltan -= n_k
                # Eliminamos los nodos nominados del dict en el cual vamos a
                # seguir buscando más nodos (ya sea para esta iteración del for
                # o para las siguientes)
                grado_to_noesenciales = elimino_nodos(grado_to_noesenciales,
                                                      k,
                                                      grado_to_noesenciales[k])
                # Actualizamos la lista de grados disponibles
                grados_disponibles = list(grado_to_noesenciales.keys())
                # Si no quedan más grados disponibles, ya está:
                if grados_disponibles == []:
                    break
                # Si todavía quedan, actualizamos el valor de k para seguir
                # buscando más nodos:
                k = obtener_mascercano(grados_disponibles, k, cercania=0)
                # cercania=0 siempre está bien porque el k viejo ya ha sido
                # eliminado de los grados disponibles
            else:
                # Con este k ya nos basta! Elijamos al azar de acá y listo.
                nominados = np.random.choice(grado_to_noesenciales[k],
                                             size=cuantos_faltan,
                                             replace=False)
                nodos_a_eliminar_temp += list(nominados)
                #Eliminamos del dict los nodos que ya marcamos para eliminar:
                grado_to_noesenciales = elimino_nodos(grado_to_noesenciales,
                                                      k,
                                                      list(nominados))
                cuantos_faltan = 0
        # Si terminó el while es porque ya marcamos suficientes nodos para
        # eliminar. Podemos pasar al siguiente for
        nodos_a_eliminar += nodos_a_eliminar_temp
    # Terminó el for. Ahora sí eliminamos todos los nodos marcados.
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
n_historias = 100 # Para las remociones aleatorias

print('--------------------------------------------------------')
print('Tamaños relativos de las componentes gigantes luego', '\n',
      'de la remoción de un subconjunto de nodos')
print('--------------------------------------------------------')
# Genero las fracciones correspondientes a remoción de esenciales
fracs_esen = np.zeros((4))
print('Eliminando nodos esenciales:')
ti = time()
for i, g in enumerate(grafos):
    fracs_esen[i] = desarme_esenciales_saco_ess(g, ess)
    print('\t', 'Grafo', i, nombres[g], ': listo')
tf = time(); print('\t', tf - ti, ' segundos')

# Genero las fracciones correspondientes a remoción de equivalentes en grado
print('Eliminando nodos no esenciales equivalentes de manera aleatoria')

# Criterio viejo (cv)
print('Criterio: "viejo" (más aleatorio)')
print('Número de historias: ', n_historias)
ti = time()
fracs_equiv_cv, sigmas_cv = np.zeros((2, 4))
for i, g in enumerate(grafos):
    resultados = np.zeros((n_historias))
    for j in range(n_historias):
        resultados[j] = efecto_nodos_equivalentes_criterioviejo(g, ess)
    fracs_equiv_cv[i] = np.average(resultados)
    sigmas_cv[i] = np.std(resultados)
    print('\t', 'Grafo', i, nombres[g], ': listo')
tf = time(); print('\t', tf - ti, ' segundos')

zscores_cv = (fracs_esen - fracs_equiv_cv) / sigmas_cv
data_cv = np.array([fracs_equiv_cv, sigmas_cv, zscores_cv])


# Criterio nuevo 
print('Criterio: "nuevo" (menos aleatorio)')
print('Número de historias: ', n_historias)
ti = time()
fracs_equiv_cn, sigmas_cn = np.zeros((2, 4))
for i, g in enumerate(grafos):
    resultados = np.zeros((n_historias))
    for j in range(n_historias):
        resultados[j] = efecto_nodos_equivalentes(g, ess)
    fracs_equiv_cn[i] = np.average(resultados)
    sigmas_cn[i] = np.std(resultados)
    print('\t', 'Grafo', i, nombres[g], ': listo')
tf = time(); print('\t', tf - ti, ' segundos')

zscores_cn = (fracs_esen - fracs_equiv_cn) / sigmas_cn
data_cn = np.array([fracs_equiv_cn, sigmas_cn, zscores_cn])

print('\n')

output_path = 'Tp2/tc02Data/datos_tabla3_100historias_2criterios.npz'
print('Guardando datos en', output_path)
np.savez(output_path,
         fracs_esen=fracs_esen,
         data_cv=data_cv,
         data_cn=data_cn)
print('Guardado con éxito!')

print('\n')

# Tabla criterio viejo
equiv_prolijo_cv = ['{:.2f} +/- {:.2f}'.format(x, y) for x, y in zip(fracs_equiv_cv, sigmas_cv)]
pd.set_option('precision', 2)
tabla3_cv = pd.DataFrame(data={'Esenciales':fracs_esen,
                            'No esenciales equivalentes:':equiv_prolijo_cv,
                            'Z-scores':zscores_cv},
                      index=['Lit', 'Lit_reg', 'APMS', 'Y2H'])

print('--------------------', '\n',
      'Tabla criterio viejo',
      '\n', '--------------------')
print(tabla3_cv)

# Tabla criterio nuevo
equiv_prolijo_cn = ['{:.2f} +/- {:.2f}'.format(x, y) for x, y in zip(fracs_equiv_cn, sigmas_cn)]
pd.set_option('precision', 2)
tabla3_cn = pd.DataFrame(data={'Esenciales':fracs_esen,
                            'No esenciales equivalentes:':equiv_prolijo_cn,
                            'Z-scores':zscores_cn},
                      index=['Lit', 'Lit_reg', 'APMS', 'Y2H'])

print('--------------------', '\n',
      'Tabla criterio nuevo',
      '\n', '--------------------')
print(tabla3_cn)

#%%
# Para importar los datos e.g. en Jupyter y graficar las tablas de forma linda

npzfile = np.load('Tp2/tc02Data/datos_tabla3_100historias_2criterios.npz')
fracs_esen = npzfile['fracs_esen']
fracs_equiv_cv, sigmas_cv, zscores_cv = npzfile['data_cv']
fracs_equiv_cn, sigmas_cn, zscores_cn = npzfile['data_cn']

equiv_prolijo_cv = ['{:.2f} +/- {:.2f}'.format(x, y) for x, y in zip(fracs_equiv_cv, sigmas_cv)]
pd.set_option('precision', 2)
tabla3_cv = pd.DataFrame(data={'Esenciales':fracs_esen,
                            'No esenciales equivalentes':equiv_prolijo_cv,
                            'Z-scores':zscores_cv},
                      index=['Lit', 'Lit_reg', 'APMS', 'Y2H'])

equiv_prolijo_cn = ['{:.2f} +/- {:.2f}'.format(x, y) for x, y in zip(fracs_equiv_cn, sigmas_cn)]
pd.set_option('precision', 2)
tabla3_cn = pd.DataFrame(data={'Esenciales':fracs_esen,
                            'No esenciales equivalentes':equiv_prolijo_cn,
                            'Z-scores':zscores_cn},
                      index=['Lit', 'Lit_reg', 'APMS', 'Y2H'])