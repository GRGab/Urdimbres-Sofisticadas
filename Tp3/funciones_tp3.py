import numpy as np
import matplotlib.pyplot as plt

import networkx as nx
from networkx import NetworkXError
from networkx.algorithms.community.community_utils import is_partition
from itertools import product
from networkx.readwrite.gml import read_gml
from histograma import histograma
import sys
import igraph as igraph

#%%

def formatear_particion(nodos, labels):
    """Dada una partición representada de manera no deseada,
    devuelve una lista de listas en la que cada sublista contiene
    los nombres de los nodos que pertenecen a una cierta comunidad.
    
    nodos : list
        Lista con los nombres de cada nodo
    labels : list
        Lista con las etiquetas de las comunidades a las que pertenece cada
        nodo.
    """
    
    output = []
    labels = np.array(labels)
    for label in set(labels):
        indices_cluster = np.where(labels == label)[0].tolist()
        cluster = [nodos[i] for i in indices_cluster]
        output.append(cluster)
    return output


def calcular_particion(nx_Graph, method="infomap", out_format='listadelistas',
                       only_labels = False):
    """
    Calcula el agrupamiento en comunidades de un grafo.
    
    In:
        nx_Graph: grafo de networkx
        method: metodo de clustering, puede ser: "infomap", "fastgreedy", "eigenvector", "louvain", "edge_betweenness","label_prop", "walktrap", ""
        
    Out:
        labels_dict: diccionario de nodo : a label al cluster al que pertenece.
    """
#    if method == "edge_betweenness":
#        nx_Graph = max(nx.connected_component_subgraphs(nx_Graph), key=len)#se queda con la componente más grande.
#        print("AVISO: restringiendo a la componente connexa más grade. De otro modo falla el algoritmo de detección de comunidades edge_betweenness.")
#    
    isdirected = nx.is_directed(nx_Graph)
    np_adj_list = nx.to_numpy_matrix(nx_Graph)
    g = igraph.Graph.Weighted_Adjacency(np_adj_list.tolist(),mode=igraph.ADJ_UPPER)
   
    if method=="infomap":
        labels = g.community_infomap(edge_weights="weight").membership
    elif method=="label_prop":
        labels = g.community_label_propagation(weights="weight").membership
    elif method=="fastgreedy":
        labels = g.community_fastgreedy(weights="weight").as_clustering().membership
    elif method=="eigenvector":
        labels = g.community_leading_eigenvector(weights="weight").membership
    elif method=="louvain":
        labels = g.community_multilevel(weights="weight").membership
    elif method=="edge_betweenness":
        labels = g.community_edge_betweenness(weights="weight", directed=isdirected).as_clustering().membership
    elif method=="walktrap":
        labels = g.community_walktrap(weights="weight").as_clustering().membership
    else:
        print('Metodo no reconocido!')
        return None
    
    if only_labels == True:
        return labels
    else:
        nodos = list(nx_Graph.nodes())
        if out_format == 'listadelistas':
            output = formatear_particion(nodos, labels)
        elif out_format == 'dict':
            output = {node:label for node,label in zip(nodos, labels)}
            
        return output

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

def comunidad_a_color(g, lista):
    """
    Funcion para asignar colores a las comunidades. Devuelve una lista con colores
    con el mismo orden que la lista de nodos, para meter en la funcion 
    nx.draw(G, node_color = comunidad_a_color(G, lista)). Además devuelve otra lista
    con los colores asociados a cada cluster.

    La lista input corresponde a la lista de listas, donde cada sublista corresponde a
    una comunidad.
    
    Input: (g, lista)   (grafo de networkx, lista de listas)
    
    Returns:
            colores_nodos : lista
                Lista con el color de cada nodo, ordenada según el orden de los nodos
                en el objeto nx.Graph
            colores_clusters : lista
                Lista con el color correspondiente a cada cluster, ordenada según el
                orden de los clusters en la lista de listas input.
    .
    .
    """
    colores_posibles = ['r', 'b', 'g', 'k', 'c', 'y', 'violet', 
                        'orange', 'indianred', 'darkgray']
    colores_random = np.random.choice(np.arange(len(colores_posibles)), size=len(lista),
                                      replace=False)
    nodos = list(g.nodes())
    colores_nodos = list(np.zeros(len(nodos)))
    colores_clusters = list(np.zeros((len(lista))))
    for i in range(len(lista)):
        index = colores_random[i]
        colores_clusters[i] = colores_posibles[index]
        for j in range(len(nodos)):
            if nodos[j] in lista[i]:
                colores_nodos[j] = colores_posibles[index]
    return colores_nodos, colores_clusters

def indices_to_nodos(graph, lista_de_indices):
    """Dada una lista de numeros enteros como indices, devuelve una lista del
    mismo tamaño reemplazando cada elemento por el nodo correspondientes a
    dicho indice."""
    lista_nueva = []
    lista_nodos = list(graph.nodes())
    for i in lista_de_indices:
        lista_nueva.append(lista_nodos[i])
    return lista_nueva


def indices_to_nodos_particion(graph, particion):
    '''Dada una lista de listas de los indicies de los nodos de una particion,
    devuelve una lista de listas del mismo tamaño reemplazando los indices
    por los nodos correspondientes.'''
    for i in range(len(particion)):
        particion[i] = indices_to_nodos(graph, particion[i])
    return particion
       
def histograma_recableo(valores, valor_obs, nombre_metodo, nombre_variable, bins=15):
    fig, ax = histograma(valores, bins=bins, density=False,
                         titulo=r'{} - Distribución de {} bajo $H_0$'
                         .format(nombre_metodo, nombre_variable),
                         xlabel=nombre_variable)
    ax.axvline(valor_obs, color='deeppink',
               label='Valor obs. = {:.2g}'.format(valor_obs))
    # Estimación del p-valor
    if valor_obs > np.mean(valores):
        pval = np.sum(np.array(valores) >= valor_obs) / len(valores)
    else:
        pval = np.sum(np.array(valores) <= valor_obs) / len(valores)
    # Si la estimación da 0, entonces decimos que p < 1/(numero de recableos)
    label = 'p = {:.2g}'.format(pval) if pval > 0 else 'p < {:2g}'.format(1/len(valores))
    ax.plot([], [], ' ', label=label)
    ax.legend(fontsize=16)
    fig.tight_layout()
    plt.show()
    return fig, ax

def histograma_recableo_multimetodo(valores, valores_obs, nombres_metodos,
                                    nombre_variable, bins=15):
    """Asume que son 8 métodos máximo."""
    n_metodos = len(nombres_metodos)
    with plt.style.context(('seaborn')):
        fig, axes = plt.subplots(2, 4, figsize=(16, 10))
        axes = np.ravel(axes)
    for i in range(n_metodos):
        histograma(valores[i], bins=bins, density=False,
                    titulo=nombres_metodos[i],
                    xlabel=nombre_variable,
                    ax=axes[i])
        axes[i].axvline(valores_obs[i], color='deeppink',
                   label='Valor obs. = {:.2g}'.format(valores_obs[i]))
        # Estimación del p-valor
        if valores_obs[i] > np.mean(valores[i]):
            pval = np.sum(np.array(valores[i]) >= valores_obs[i]) / len(valores[i])
        else:
            pval = np.sum(np.array(valores[i]) <= valores_obs[i]) / len(valores[i])
        # Si la estimación da 0, entonces decimos que p < 1/(numero de recableos)
        label = 'p = {:.2g}'.format(pval) if pval > 0 else 'p < {:2g}'.format(1/len(valores[i]))
        axes[i].plot([], [], ' ', label=label)
        axes[i].legend(fontsize=12)
    for i in range(n_metodos, 8):
        axes[i].axis('off')
    plt.show()
    return fig, axes

#%%
if __name__ == '__main__':
    import matplotlib.pyplot as plt
    plt.ion()
    import time
    sys.path.append('./Tp3/')
    from silhouettes import silhouettes

    G = nx.balanced_tree(h=3,r=2)
#    nx.draw(G,with_labels=True)
#    plt.show()
    print('Prueba con infomap')
    particion = calcular_particion(G, method='infomap')
    modularidad = calcular_modularidad(G, particion)
    print('La modularidad es', modularidad)
    colores, _ = comunidad_a_color(G, particion)

