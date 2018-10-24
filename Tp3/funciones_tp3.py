# -*- coding: utf-8 -*-
"""
Created on Fri Oct 19 19:13:06 2018

@author: Gabo
"""

import matplotlib.pyplot as plt
import numpy as np

import time

import networkx as nx
from networkx.algorithms.community.centrality import girvan_newman
from networkx import NetworkXError
from networkx.utils import not_implemented_for
from networkx.algorithms.community.community_utils import is_partition
from networkx.readwrite.gml import read_gml


import igraph as igraph
import community # instalar como python-Louvain


#%%

def calculate_partition(np_adj_list, method="infomap"):
    t0=time.time()
    if method in ['infomap', 'fastgreedy']:
        g = igraph.Graph.Weighted_Adjacency(np_adj_list.tolist(),mode=igraph.ADJ_UPPER)
        if method=="infomap":
            labels = g.community_infomap(edge_weights="weight").membership
    #    labels = g.community_label_propagation(weights="weight").membership
        if method=="fastgreedy":
            labels = g.community_fastgreedy(edge_weights="weight").membership
    if method == 'edge_bet':
       labels = girvan_newman(g)
    print("Duración: {}s".format(time.time()-t0))
    return labels


__all__ = ['coverage', 'modularity', 'performance'] # ??????????????


class NotAPartition(NetworkXError):
    """Raised if a given collection is not a partition.

    """

    def __init__(self, G, collection):
        msg = '{} is not a valid partition of the graph {}'
        msg = msg.format(G, collection)
        super(NotAPartition, self).__init__(msg)


def modularity(G, communities, weight='weight'):
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
def atribuir_comunidades(criterio):
    """Toma una lista con las etiquetas de las comunidades a las que pertenece
    cada nodo, y devuelve una lista de listas en la que cada sublista contiene
    a todos los nodos de una única comunidad."""
    lista = []
    for i in np.unique(criterio):
        lista_por_comunidades = np.where(np.array(criterio) == i)[0].tolist()
        if len(lista_por_comunidades) > 0:
            lista.append(lista_por_comunidades)
    return lista

#%%
""" Hacer particion con louvain"""
partition = community.best_partition(G)
#%%
""" fast-greedy sin igraph"""
import networkx as nx
import matplotlib.pyplot as plt
G = nx.balanced_tree(h=3,r=2)
nx.draw(G,with_labels=True)
plt.show()
#%%
comus = nx.algorithms.community.greedy_modularity_communities(G, weight=None)
list(comus)
#%%
"""Girvan-newman"""
com = nx.algorithms.community.centrality.girvan_newman(G)
com = list(com)
com = [list([list(conjunto) for conjunto in tupla]) for tupla in com]
