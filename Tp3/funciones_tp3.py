# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 19:13:06 2018

@author: Gabo
"""

import numpy as np

import networkx as nx
from networkx import NetworkXError
from networkx.algorithms.community.community_utils import is_partition
from networkx.readwrite.gml import read_gml


import igraph as igraph

#%%
def formatear_particion(particion):
    """Dada una partición representada de alguna manera no deseada,
    devuelve una lista de listas en la que cada sublista contiene
    el número de cada uno de los nodos que pertenece a una cierta comunidad.
    
    'particion' puede ser una lista con las etiquetas de las comunidades a las
    que pertenece cada nodo, o bien un diccionario que a cada nodo le asigna
    el número de su cluster.
    
    Todo esto asume que los nodos están bien ordenados siempre.
    
    PENDIENTE
    ---------
    - Reemplazar los mensajes de error por errores posta."""
    
    output = []
    if isinstance(particion, list):
        if isinstance(particion[0], list):
            print('La partición ya pareciera estar en el formato deseado.')
            return
        else:
            for i in set(particion):
                lista_por_comunidades = np.where(np.array(particion) == i)[0].tolist()
                output.append(lista_por_comunidades)
    elif isinstance(particion, dict):
        for i in set(particion.values()):
            lista_por_comunidades = []
            for j in particion.keys():
                if particion[j] == i:
                    lista_por_comunidades.append(j)
            output.append(lista_por_comunidades)
    else:
        print('No se reconoce el formato de entrada.')
        return
    return output


def calcular_particion(nx_Graph, method="infomap", out_format='listadelistas'):
    """
    Calcula el agrupamiento en comunidades de un grafo.
    
    In:
        nx_Graph: grafo de networkx
        method: metodo de clustering, puede ser: "infomap", "fastgreedy", "eigenvector", "louvain", "edge_betweenness","label_prop", "walktrap", ""
        
    Out:
        labels_dict: diccionario de nodo : a label al cluster al que pertenece.
    """
    if method == "edge_betweenness":
        nx_Graph = max(nx.connected_component_subgraphs(nx_Graph), key=len)#se queda con la componente más grande.
        print("AVISO: restringiendo a la componente connexa más grade. De otro modo falla el algoritmo de detección de comunidades edge_betweenness.")
    
    isdirected = nx.is_directed(nx_Graph)
    np_adj_list = nx.to_numpy_matrix(nx_Graph)
    g = igraph.Graph.Weighted_Adjacency(np_adj_list.tolist(),mode=igraph.ADJ_UPPER)
   
    if method=="infomap":
        labels = g.community_infomap(edge_weights="weight").membership
    if method=="label_prop":
        labels = g.community_label_propagation(weights="weight").membership
    if method=="fastgreedy":
        labels = g.community_fastgreedy(weights="weight").as_clustering().membership
    if method=="eigenvector":
        labels = g.community_leading_eigenvector(weights="weight").membership
    if method=="louvain":
        labels = g.community_multilevel(weights="weight").membership
    if method=="edge_betweenness":
        labels = g.community_edge_betweenness(weights="weight", directed=isdirected).as_clustering().membership
    if method=="walktrap":
        labels = g.community_walktrap(weights="weight").as_clustering().membership
    
    if out_format == 'listadelistas':
        output = formatear_particion(labels)
    elif out_format == 'dict':
        output = {node:label for node,label in zip(nx_Graph.nodes(), labels)}
    return output


#%%
class NotAPartition(NetworkXError):
    """Raised if a given collection is not a partition.

    """

    def __init__(self, G, collection):
        msg = '{} is not a valid partition of the graph {}'
        msg = msg.format(G, collection)
        super(NotAPartition, self).__init__(msg)


def calcular_modularidad(G, communities, weight='weight'):
    r"""Returns the modularity of the given partition of the graph.

    Modularity is defined in [1]_ as

    .. math::

        Q = \frac{1}{2m} \sum_{ij} \left( A_{ij} - \frac{k_ik_j}{2m}\right)
            \delta(c_i,c_j)

    where $m$ is the number of edges, $A$ is the adjacency matrix of
    `G`, $k_i$ is the degree of $i$ and $\delta(c_i, c_j)$
    is 1 if $i$ and $j$ are in the same community and 0 otherwise.

    Parameters
    ----------
    G : NetworkX Graph

    communities : list
        List of sets of nodes of `G` representing a partition of the
        nodes.

    Returns
    -------
    Q : float
        The modularity of the paritition.

    Raises
    ------
    NotAPartition
        If `communities` is not a partition of the nodes of `G`.

    Examples
    --------
    >>> G = nx.barbell_graph(3, 0)
    >>> nx.algorithms.community.modularity(G, [{0, 1, 2}, {3, 4, 5}])
    0.35714285714285704

    References
    ----------
    .. [1] M. E. J. Newman *Networks: An Introduction*, page 224.
       Oxford University Press, 2011.

    """
    if not is_partition(G, communities):
        raise NotAPartition(G, communities)

    multigraph = G.is_multigraph()
    directed = G.is_directed()
    m = G.size(weight=weight)
    if directed:
        out_degree = dict(G.out_degree(weight=weight))
        in_degree = dict(G.in_degree(weight=weight))
        norm = 1 / m
    else:
        out_degree = dict(G.degree(weight=weight))
        in_degree = out_degree
        norm = 1 / (2 * m)

    def val(u, v):
        try:
            if multigraph:
                w = sum(d.get(weight, 1) for k, d in G[u][v].items())
            else:
                w = G[u][v].get(weight, 1)
        except KeyError:
            w = 0
        # Double count self-loops if the graph is undirected.
        if u == v and not directed:
            w *= 2
        return w - in_degree[u] * out_degree[v] * norm

    Q = sum(val(u, v) for c in communities for u, v in product(c, repeat=2))
    return Q * norm

#%%
if __name__ == '__main__':
    import matplotlib.pyplot as plt
    import time

    G = nx.balanced_tree(h=3,r=2)
    nx.draw(G,with_labels=True)
    plt.show()
    
    particion = calcular_particion(G, method='infomap')
    modularidad = calcular_modularidad(G, particion)
    print('La modularidad es', modularidad)
