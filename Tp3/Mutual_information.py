import numpy as np
import networkx as nx
from networkx.readwrite.gml import read_gml
from lectura import ldata
import sys
sys.path.append('./Tp3/')
from funciones_tp3 import calcular_particion, comunidad_a_color
import matplotlib.pyplot as plt
import igraph as igraph
#%%
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

def I_M(g, method_1, method_2, genero = False):
    colores_posibles = ['r', 'b', 'g', 'k', 'c', 'y', 'violet','orange', 
                        'indianred', 'darkgray']
    
    if genero:
        nodes_1 = particionar_por_genero(g)
        colors_1 = comunidad_a_color(g, nodes_1)  
    else:
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
if __name__ == '__main__':
    dolph = read_gml('Tp3/dolphins.gml')    
    lista = ["infomap","label_prop", "fastgreedy", "eigenvector", "louvain"
            , "edge_betweenness", "walktrap"]
    
    metodos = ["fastgreedy", "eigenvector","edge_betweenness", "louvain", "walktrap", "infomap", "label_prop"]
    
    
    tabla = np.zeros([len(metodos), len(metodos)])
    for i in range(len(metodos)):
        for j in range(len(metodos)):
            tabla[i, j] = round(I_M(dolph, metodos[i], metodos[j]), 4)
    
    import copy
    tabla2 = copy.deepcopy(tabla)
    
    tabla2 = (tabla + np.transpose(tabla)) / 2
    
    import pandas as pd
    cuadro = pd.DataFrame(tabla2, columns = metodos, index = metodos)
#%%
    genders = dict(ldata('Tp3/dolphinsGender.txt'))
    
    # Agrego los sexos a los dicts de cada delfín
    for nodo, dict_nodo in dict(dolph.nodes).items():
        dict_nodo['gender'] = genders[nodo] # agrego el sexo del delfín a su dict
    #    print('Key = {}, Value = {}'.format(nodo, dict_nodo)) # para chequear que anda
    
    tabla = np.zeros(len(metodos))
    for i in range(len(metodos)):
        tabla[i] = I_M(dolph, 'asd', metodos[i], genero = True)
    np.save('tabla_info_mutua_generos', tabla)
    import pandas as pd
    cuadro = pd.DataFrame(tabla, columns = ['Genero'], index = metodos)
    
    