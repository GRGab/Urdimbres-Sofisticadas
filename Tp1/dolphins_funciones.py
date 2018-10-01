# -*- coding: utf-8 -*-
"""
Created on Sun Sep  2 21:13:35 2018

@author: Gabo
"""

import networkx as nx
import matplotlib.pyplot as plt
from collections import Counter
import numpy as np

################## Para graficar 
def genero_a_color(gender):
    if gender=='m':
        return 'red'
    elif gender=='f':
        return 'dodgerblue'
    else:
        return 'green'
    
def particionar_por_genero(G, orden={'f':0, 'NA':1, 'm':2}):
    particiones = [[], [], []]
    for key in dict(G.nodes).keys():
        gender = G.nodes[key]['gender']
        if gender=='f':
            particiones[orden[gender]].append(key)
        elif gender=='m':
            particiones[orden[gender]].append(key)
        else:
            particiones[orden[gender]].append(key)
    return particiones

def crear_leyenda(ax):
    ax.plot([],[], 'o', color='dodgerblue', label='Hembras')
    ax.plot([],[], 'o', color='red', label='Machos')
    ax.plot([],[], 'o', color='black', label='NA')
    ax.legend()
    
################## Para cuantificar homofilia

def contar_clases(g, atributo, valores_posibles):
    ns = []
    for valor in valores_posibles:
        n = len([n for n, attrdict in dict(g.nodes).items() if attrdict[atributo]==valor])
        ns.append(n)
        print('Hay {} nodos con {}={}'.format(n, atributo, valor))
    return ns

def contar_enlaces_internos(g, atributo, valor):
    """Cuenta los enlaces internos en el grupo de nodos que tienen
    atributo=valor. g debe ser objeto Graph de NetworkX.
    Ejemplo con delfines: atributo='gender', valor puede ser 'f', 'm' o 'NA'
    """
    nodos = [n for n, attrdict in dict(g.nodes).items() if attrdict[atributo]==valor]
    grupo = nx.subgraph(g, nodos)
    return grupo.size()

def contar_enlaces_entre_grupos(g, atributo):
    """Cuenta los enlaces que conectan grupos con valores distintos de
    atributo. g debe ser objeto Graph de NetworkX.
    Ejemplo con delfines: atributo='gender'.
    """
    n = 0
    for edge in g.edges():
        a, b = edge[0], edge[1]
        if g.nodes()[a][atributo] != g.nodes()[b][atributo]:
            n = n + 1
    return n

################## Calculo del p-value
def p_value(datos, valor_real=52):
    enlace = list(Counter(datos).keys())
    cuantos = list(Counter(datos).values())
    indices = np.where(np.array(enlace)<= valor_real)[0] 
    a = [cuantos[i] for i in indices]
    return 2 * len(a)/len(datos)


################## Para desarmar la red
    
#Esta funcion agarra un grafo, y lo va desarmando sacando una n cantidad de nodos
#con el criterio de que esos nodos sean los de coeficiente de clustering mas chico
#de la red.     La funcion tiene input el grafo (g), y el porcentaje de nodos que 
#queres sacarle a la red (porcentaje_nodos_sacados). El porcentaje se escribe
#como 0.1 para sacar el 10% de los nodos.

def desarme(g, porcentaje_nodos_sacados):
    n= int(len(g.nodes.keys()) * porcentaje_nodos_sacados)
    G = g.copy()
    nodes = []
    clustering = []
    for a in G.nodes.keys():
        if G.degree(a)>1:
            nodes.append(a)
            clustering.append(nx.clustering(G)[a])
    for i in range(n):
        j = clustering.index(min(clustering))
        G.remove_node(nodes[j])
        clustering.pop(j)
        nodes.pop(j)
    return G

def cociente(g):
    '''Agarra una red G y va sacando una cierta proporcion de nodos hasta que
    el cociente entre la cant de nodos de la componente gigante y la segunda 
    comoponente sea cercano a 1.'''
    rango_porcentajes = np.arange(0, 0.85, 1/g.order())
    d = []
    for i in rango_porcentajes:
        g_desarmado = desarme(g, i)
        lengths = [len(d) for d in sorted(nx.connected_components(g_desarmado), key=len, reverse=True)]
        if len(lengths)>=2:
            d.append(lengths[0]/lengths[1])
    d2 = np.asarray(d)
    d2 = abs(d2-1)
    d2 = list(d2)
    index = d2.index(min(d2))
    porcentaje = rango_porcentajes[index]
    return porcentaje, min(d)
