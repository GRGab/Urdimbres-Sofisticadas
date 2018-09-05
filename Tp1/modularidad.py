import networkx as nx
import matplotlib.pyplot as plt
 
def modularidad(g, atributo):
    A = nx.adjacency_matrix(g)
    B = 0
    m = len(g.edges())
    k = []
    name = []
    for key, value in g.degree():
        k.append(value)
        name.append(key)
    for i in range(len(k)):
        for j in range(len(k)):
            if g.node[name[i]][atributo] == g.node[name[j]][atributo]:
                B += A[i, j] - k[i]*k[j]/(2*m)
    return float(B)/(2*m)

if __name__ == '__main__':
    g = nx.gnp_random_graph(10, 0.3)
    for node, attrDict in dict(g.nodes()).items():
        attrDict['gender'] = 'a' if (node % 2 == 0) else 'b'
    
    colores = [('blue' if gender=='a' else 'red') for gender in nx.get_node_attributes(g, "gender").values()]
    plt.figure()
    nx.draw(g, node_color=colores)
    print(modularidad(g, 'gender'))