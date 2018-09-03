# -*- coding: utf-8 -*-
"""
Created on Sun Sep  2 21:13:35 2018

@author: Gabo
"""

import networkx as nx
import matplotlib.pyplot as plt

################## Para graficar 
def genero_a_color(gender):
    if gender=='m':
        return 'red'
    elif gender=='f':
        return 'dodgerblue'
    else:
        return 'black'
    
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