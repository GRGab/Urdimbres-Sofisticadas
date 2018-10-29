import numpy as np
import networkx as nx
from networkx.readwrite.gml import read_gml
import sys
sys.path.append('./Tp3/')
from funciones_tp3 import calcular_particion, comunidad_a_color
import matplotlib.pyplot as plt
import igraph as igraph
#%%
dolph = read_gml('Tp3/dolphins.gml')    
lista = ["infomap","label_prop", "fastgreedy", "eigenvector", "louvain"
        , "edge_betweenness", "walktrap"]

metodos = ["infomap", "fastgreedy", "eigenvector", "louvain", "edge_betweenness","label_prop", "walktrap"]


def I_M(g, method_1, method_2):
    colores_posibles = ['r', 'b', 'g', 'k', 'c', 'y', 'violet','orange', 
                        'indianred', 'darkgray']
    
    nodes_1 = calcular_particion(g, method = method_1)
    colors_1 = comunidad_a_color(g, nodes_1)
    
    nodes_2 = calcular_particion(g, method = method_2)
    colors_2 = comunidad_a_color(g, nodes_2)
    
    matrix = np.zeros([len(colores_posibles), len(colores_posibles)])
    for i in range(len(g.nodes())):
        index_1 = colores_posibles.index(colors_1[0][i])
        index_2 = colores_posibles.index(colors_2[0][i])
        matrix[index_1, index_2] += 1
    matrix = matrix / len(g.nodes())
    
    colors_1 = [colores_posibles.index(colors_1[0][i]) for i in range(len(colors_1[0]))]
    colors_2 = [colores_posibles.index(colors_2[0][i]) for i in range(len(colors_2[0]))]
    m_1, bins = np.histogram(colors_1, bins = np.arange(len(colores_posibles)+1))
    m_2, bins = np.histogram(colors_2, bins = np.arange(len(colores_posibles)+1))
    m_1 = m_1 / len(g.nodes())
    m_2 = m_2 / len(g.nodes())
    
    I = 0
    for i in range(len(colores_posibles)):
        for j in range(len(colores_posibles)):
            if matrix[i, j] != 0 and m_1[i] != 0 and m_2[j] != 0:
                I += matrix[i, j] * np.log2(matrix[i, j]/(m_1[i] * m_2[j]))
                
    H_1 = - np.sum([i * np.log2(i) for i in m_1 if i != 0])
    H_2 = - np.sum([i * np.log2(i) for i in m_2 if i != 0]) 
    
    I_M = I / (0.5 * (H_1 + H_2))
    
    return I_M
#%%

tabla = np.zeros([len(metodos), len(metodos)])
for i in range(len(metodos)):
    for j in range(len(metodos)):
        tabla[i, j] = round(I_M(dolph, metodos[i], metodos[j]), 5)

