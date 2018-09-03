import networkx as nx

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
        nodes.append(a)
        clustering.append(nx.clustering(G)[a])
    for i in range(n):
        j = clustering.index(min(clustering))
        G.remove_node(nodes[j])
        clustering.pop(j)
        nodes.pop(j)
    return G


#Esto de aca abajo es una forma de ver como se destruye la red dolphin para distintos
#valores de porcentaje_nodos_sacados.
#plt.figure(1)
#nx.draw(dolph, node_size=30)
#
#for i in np.arange(0, 0.8, 0.1):
#    grafo = desarme(dolph, i)
#    plt.figure(i*10 + 1)
#    nx.draw(grafo, node_size=30)
# 