import networkx as nx
import numpy as np
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

print(cociente(dolph))


#Esto de aca abajo es una forma de ver como se destruye la red dolphin para distintos
#valores de porcentaje_nodos_sacados.
#plt.figure(1)
#nx.draw(dolph, node_size=30)
##
#for i in np.arange(0, 0.8, 0.1):
#    grafo = desarme(dolph, i)
#    plt.figure(i*10 + 1)
#    nx.draw(grafo, node_size=30)
 