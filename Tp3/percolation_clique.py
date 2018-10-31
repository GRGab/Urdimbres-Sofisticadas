import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from networkx.readwrite.gml import read_gml
from collections import defaultdict

dolph = read_gml('dolphins.gml')    

#genders = dict(ldata('Tp3/dolphinsGender.txt'))
#
## Agrego los sexos a los dicts de cada delfín
#for nodo, dict_nodo in dict(dolph.nodes).items():
#    dict_nodo['gender'] = genders[nodo]

def k_clique_communities(G, k, cliques=None):
    """Find k-clique communities in graph using the percolation method.

    A k-clique community is the union of all cliques of size k that
    can be reached through adjacent (sharing k-1 nodes) k-cliques.

    Parameters
    ----------
    G : NetworkX graph

    k : int
       Size of smallest clique

    cliques: list or generator       
       Precomputed cliques (use networkx.find_cliques(G))

    Returns
    -------
    Yields sets of nodes, one for each k-clique community.

    Examples
    --------
    >>> G = nx.complete_graph(5)
    >>> K5 = nx.convert_node_labels_to_integers(G,first_label=2)
    >>> G.add_edges_from(K5.edges())
    >>> c = list(nx.k_clique_communities(G, 4))
    >>> list(c[0])
    [0, 1, 2, 3, 4, 5, 6]
    >>> list(nx.k_clique_communities(G, 6))
    []

    References
    ----------
    .. [1] Gergely Palla, Imre Derényi, Illés Farkas1, and Tamás Vicsek,
       Uncovering the overlapping community structure of complex networks 
       in nature and society Nature 435, 814-818, 2005,
       doi:10.1038/nature03607
    """
    if k < 2:
        raise nx.NetworkXError("k=%d, k must be greater than 1."%k)
    if cliques is None:
        cliques = nx.find_cliques(G)
    cliques = [frozenset(c) for c in cliques if len(c) >= k]

    # First index which nodes are in which cliques
    membership_dict = defaultdict(list)
    for clique in cliques:
        for node in clique:
            membership_dict[node].append(clique)

    # For each clique, see which adjacent cliques percolate
    perc_graph = nx.Graph()
    perc_graph.add_nodes_from(cliques)
    for clique in cliques:
        for adj_clique in _get_adjacent_cliques(clique, membership_dict):
            if len(clique.intersection(adj_clique)) >= (k - 1):
                perc_graph.add_edge(clique, adj_clique)

    # Connected components of clique graph with perc edges
    # are the percolated cliques
    for component in nx.connected_components(perc_graph):
        yield(frozenset.union(*component))

def _get_adjacent_cliques(clique, membership_dict):
    adjacent_cliques = set()
    for n in clique:
        for adj_clique in membership_dict[n]:
            if clique != adj_clique:
                adjacent_cliques.add(adj_clique)
    return adjacent_cliques

#%%
colors = []

for j in [2, 3, 4, 5, 6, 7, 8]:
    c = list(k_clique_communities(dolph, j))
    for i in range(len(c)):
        c[i] = list(c[i])
        
    colores = comunidad_a_color(dolph, c)
    #le pongo color aqua a los que no tienen comunidad.
    for i in range(len(colores[0])):
        if colores[0][i] == 0.0:
            colores[0][i] = 'aqua'
            
    colors.append(colores[0])



#[ax1, ax2, ax3, ax4], [ax5, ax6, ax7, ax8]

ns = 35
#nx.draw(dolph, node_color = colors[0])
pos = nx.spring_layout(dolph)

fig, ([ax1, ax2, ax3], [ax4, ax5, ax6]) = plt.subplots(2, 3)

plt.sca(ax1)
ax1.set_title('k = 2')
nx.draw(dolph, node_size = ns, node_color=colors[0], pos = pos)

plt.sca(ax2)
ax2.set_title('k = 3')
nx.draw(dolph, node_size = ns, node_color=colors[1], pos = pos)

plt.sca(ax3)
ax3.set_title('k = 4')
nx.draw(dolph, node_size = ns, node_color=colors[2], pos = pos)

plt.sca(ax4)
ax4.set_title('k = 5')
nx.draw(dolph, node_size = ns, node_color=colors[3], pos = pos)

plt.sca(ax5)
ax5.set_title('k = 6')
nx.draw(dolph, node_size = ns, node_color=colors[4], pos = pos)

plt.sca(ax6)
ax6.set_title('k = 7')
nx.draw(dolph, node_size = ns, node_color=colors[5], pos = pos)

