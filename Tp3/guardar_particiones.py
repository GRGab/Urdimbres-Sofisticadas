import numpy as np
import matplotlib.pyplot as plt

import networkx as nx
from networkx import NetworkXError
from networkx.algorithms.community.community_utils import is_partition
from itertools import product
from networkx.readwrite.gml import read_gml
from histograma import histograma
import sys
sys.path.append('./Tp3/')
import igraph as igraph
from silhouettes import silhouettes
from funciones_tp3 import * 
#%%
def guardar_particiones(graph_original, N_swaps,
                        Numero_de_recableos, lista_de_metodos,
                        guardar_grafos = False, 
                        output_path = None,
                        silencioso = False):
    
    ''' Toma el grafo orginal y realiza N recableos de la red. Para cada 
    recableo, calcula las particion de la red segun un metodo de particion
    (por ejemplo, infomap) utilizando la funcion calcular_particiones. Si 
    guardar_grafo = True, la funcion guarda la informacion en un array de N*M,
    donde N es la cant de recableos y M es la cant de metodos para particionar.
    Cada elemento del array es una lista de listas con los clusters de la red.
    Ademas, guarda un array donde el i-esimo elemento es la i -esima red
    de networkx recableada. 
    Finalmente, para el grafo original, devuelve una lista con M elementos, 
    donde cada uno es una lista de listas para cada particion.
    
    Si guardar_grafos = False, las modularidades se guardan en un array de M*N.
    Ademas, e guarda un array de M componentes para todos los metodo aplicados 
    a la red original.
    Para siluette, tambien devuelve una lista, donde cada elemento representa
    la lista con los siluette de cada nodo para una particion. Es decir, 
    la longitud de la lista total es N, siendo N la cant de recableos. 
    Finalmente, tambien devuelve la lista con los siluettes del grafo original.

    Guarda
    ------
    mod_original : lista
        Cada elemento es un valor de modularidad según algún método, en el orden del
        input lista_de_metodos.
    mod_rewire : ndarray
        Array de dimensión 2 y shape (n_metodos, n_recableos). Cada elemento es un
        valor de modularidad.
    sil_original : listas anidadas (3 capas)
        Cada elemento es una lista de listas correspondiente a la partición obtenida
        mediante algún método, en el orden del input lista_de_metodos. En cada partición,
        en el lugar de cada nodo se encuentra su valor de silhouette.
    sil_rewire : lista anidadas (4 capas)
        Estructura de listas anidadas en 4 capas. De afuera hacia adentro, las capas son:
        método, recableo, partición y nodo.
        

    '''
    if guardar_grafos == False:
        ### RED ORIGINAL
        n_metodos = len(lista_de_metodos)
        mod_original = np.zeros((n_metodos))
        sil_original = []

        for i, metodo in enumerate(lista_de_metodos):
            particion = calcular_particion(graph_original, method=metodo)
            mod = calcular_modularidad(graph_original, particion)
            mod_original[i] = mod
            sil_actual = silhouettes(graph_original, particion, silencioso=silencioso)
            sil_original.append(sil_actual)

        ### RECABLEO
        # print('Comienza el recableo') # debugging
        mod_rewire = np.zeros((n_metodos, Numero_de_recableos))
        sil_rewire = [[] for _ in range(n_metodos)]
        G = graph_original.copy()
        for i in range(Numero_de_recableos): 
            # Movemos de lugar los enlaces una buena cantidad de veces
            G = nx.double_edge_swap(G, nswap=N_swaps, max_tries=N_swaps * 1.5)
            # Nos fijamos si el grafo final es conexo. Si no lo es,
            # seguimos swapeando hasta que lo sea
            while not nx.is_connected(G):
                G = nx.double_edge_swap(G, nswap=1, max_tries=10)
            # Ahora sí, calculamos las particiones y sus respectivos
            # modularidad y silhouettes
            for j, metodo in enumerate(lista_de_metodos):
                particiones_rewire = calcular_particion(G, method=metodo)
                mod_rewire[j, i] = calcular_modularidad(G, particiones_rewire)                 
                sil_rewire[j].append(silhouettes(G, particiones_rewire,
                                                 silencioso=silencioso))
            
        #Guardamos
        if output_path == None:
            output_path = 'Tp3/tc03Data/Ej_b_particiones_numeros.npz'
        np.savez(output_path, mod_original=mod_original, 
                              mod_rewire = mod_rewire,
                              sil_original = sil_original,
                              sil_rewire = sil_rewire)
        
    if guardar_grafos == True: 
        ### ATENCIÓN: CÓDIGO DESACTUALIZADO
        salida_grafo_original = [] 
        for metodo in lista_de_metodos:
            lista_nodos_original = calcular_particion(graph_original,
                                                      method=metodo)
            salida_grafo_original.append(lista_nodos_original)
        G = graph_original.copy()    
        salida = []
        grafos_rewire = []
        for metodo in lista_de_metodos:
            salida_por_rewire = []
            for i in range(Numero_de_recableos): 
    
                g_rewire = nx.double_edge_swap(G, nswap=N_swaps,
                                               max_tries=N_swaps * 1.5)
                lista_nodos = calcular_particion(g_rewire, method = metodo)
                salida_por_rewire.append(lista_nodos)
                grafos_rewire.append(g_rewire)
            salida.append(salida_por_rewire)
        if output_path == None:
            output_path = 'Tp3/tc03Data/Ej_b_con_grafos.npz'
        np.savez(output_path, salida = salida, 
                 salida_grafo_original = salida_grafo_original,
                 grafos_rewire = grafos_rewire) 
#%%
# if __name__ == '__main__':
#     import time
#     sys.path.append('./Tp3/')
#     G = nx.balanced_tree(h=5,r=2)
#     lista_de_metodos = ["infomap", "label_prop"] 
# #    guardar_particiones(G, 200, 10 , lista_de_metodos,
# #                        output_path = 'Tp3/tc03Data/pruebita_Ej_b.npz')
#     dolph = read_gml('Tp3/dolphins.gml')
#     guardar_particiones(dolph, 100, 10, lista_de_metodos,
#                         output_path='lala')
#     #%%Importamos
#     npzfile = np.load('Tp3/tc03Data/pruebita_Ej_b.npz')
#     rewire = npzfile['mod_rewire']
#     original = npzfile['mod_original']