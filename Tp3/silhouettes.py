import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import networkx as nx
from networkx.algorithms.community.community_utils import is_partition

from networkx.readwrite.gml import read_gml


import sys
sys.path.append('./Tp3')
from funciones_tp3 import (calcular_particion, NotAPartition, indices_to_nodos_particion,
                           comunidad_a_color)

def crear_nodos_to_indices(particion):
    """Crea un diccionario que a cada nombre de nodo le asigna los índices
    (m, n) tales que particion[m][n] == nombre.
    
    Input
    -----
    particion : list
        lista de listas. Cada sublista es un cluster y sus elementos son los
        nombres de los nodos que pertenecen a dicho cluster.
    Output
    -----
    out : dict
        keys son los nombres de los nodos, values son sus índices en la partición
    """
    dicc = {}
    for m in range(len(particion)):
        for n in range(len(particion[m])):
            dicc[particion[m][n]] = (m, n)
    return dicc

def silhouettes(G, particion):
    """
    Calcula el valor de silhouette para cada nodo del grafo 'G' dada una
    partición 'particion' como lista de listas. Dicho valor está dado por
    
    s(i) = (b(i) - a(i)) / max(a(i), b(i))
    
    donde a(i) es la distancia media a todos los nodos del mismo cluster que i
    y b(i) es la mínima de las distancias medias a los distintos clusters a los
    cuales no pertenece i. Para mayor claridad, sea c_i el cluster al que
    pertenece i, y sea Q = particion - c_i el conjunto de los clusters a los cuales
    no pertenece i. Entonces se define
    
    b(i) = min{promedio{d(i,j) : j in cluster} : cluster in Q}
    
    b(i) también se suele llamar "distancia media al cluster más cercano".

    Input
    -----
    G : nx.Graph
    particion : list
        lista de listas. Cada sublista es un cluster y sus elementos son los
        nombres de los nodos que pertenecen a dicho cluster.
    Output
    ------
    output : list
        lista de listas. Cada sublista es un cluster y sus elementos son los
        valores de silhouette para cada nodo, preservando el orden del input.
    """
    if not is_partition(G, particion):
        raise NotAPartition(G, particion)

    ds = list(nx.all_pairs_shortest_path_length(G))
    d = lambda i,j: ds[i][1][j]
    # ds[i][1][j] es la distancia (longitud del camino más corto)
    # entre i y j

    n = G.order()
    nc = len(particion)
    # Creamos lista de lista con iguales longitudes que 'particion'
    s_values = [[[] for n in range(len(particion[m]))] for m in range(nc)]
    # La iremos rellenando con valores de silhouette.
    nodos_to_indices = crear_nodos_to_indices(particion)
    # Recorremos los nodos en el ordenamiento global correspondiente
    # a la función distancia 'd'
    for i, nodo in enumerate(G.nodes()):
        m, n = nodos_to_indices[nodo]
        cluster_actual = particion[m]
        otros_clusters = (particion[l] for l in range(nc) if l != m)
        a = np.average([d(i,j) for j in cluster_actual])
        
        dists_interclusters = [np.average([d(i,j) for j in cluster]) \
                                                  for cluster in otros_clusters]
        b = min(dists_interclusters)
        s_values[m][n] = (b - a) / max(a, b)
    return s_values

def graficar_silhouettes(sil_vals, colores=None, ax=None, titulo=None):
    n_clusters = len(sil_vals)
    sils_flattened = [s for cluster in sil_vals for s in cluster]
    sil_medio = np.mean(sils_flattened)
    n_nodos = len(sils_flattened)
    # The silhouette coefficient can range from -1, 1
    
    if ax is None:
        with plt.style.context(('seaborn')):
            fig, ax = plt.subplots()
    else:
        fig = ax.get_figure()
    
    xmin = min(0, min(sils_flattened) - 0.01)
    xmax = max(sils_flattened) + 0.01
    ax.set_xlim([xmin, xmax])

    # The (n_clusters+1)*sep is for inserting blank space between silhouette
    # plots of individual clusters, to demarcate them clearly.
    sep = 0
    ax.set_ylim([0, n_nodos + (n_clusters + 1) * sep])
    ax.set_yticks([])  # Clear the yaxis labels / ticks
    
    y_lower = sep
    for i in range(n_clusters):
        # Aggregate the silhouette scores for samples belonging to
        # cluster i, and sort them
        ith_cluster_silhouette_values = sorted(sil_vals[i])

        size_cluster_i = len(ith_cluster_silhouette_values)
        y_upper = y_lower + size_cluster_i

        if colores is None:
            color = cm.nipy_spectral(float(i) / n_clusters)
        else:
            color = colores[i]
        ax.fill_betweenx(np.arange(y_lower, y_upper),
                          0, ith_cluster_silhouette_values,
                          facecolor=color, edgecolor=color, alpha=0.7)

        # Label the silhouette plots with their cluster numbers at the middle
        ax.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))

        # Compute the new y_lower for next plot
        y_lower = y_upper + sep

    if titulo is not None:
        ax.set_title(titulo)
    ax.set_xlabel("Coeficiente de Silhouette")
    ax.set_ylabel("Comunidad")

    # The vertical line for average silhouette score of all the values
    ax.axvline(x=sil_medio, color="red", linestyle="--")

    if ax is None:
        fig.tight_layout()
    fig.show()
    
#%%
if __name__ == '__main__':
    plt.ion()
    # G = nx.balanced_tree(h=3,r=2)
    # particion = calcular_particion(G, method='infomap')
    # sil = silhouettes(G, particion)
    #%%
    dolph = read_gml('Tp3/dolphins.gml')
    npzfile = np.load('Tp3/tc03Data/Ej_b_particiones.npz')
    rewire = npzfile['salida']
    original = npzfile['salida_grafo_original']
    particion = original[0]
    sil = silhouettes(dolph, particion)

    colores_nodos, colores_clusters = comunidad_a_color(dolph, particion)
    with plt.style.context(('seaborn')):
        fig, (ax1, ax2) = plt.subplots(1, 2)#, figsize=(10, 8))
    nx.draw(dolph, node_color=colores_nodos, ax=ax1, node_size=50)
    graficar_silhouettes(sil, ax=ax2, colores=colores_clusters)