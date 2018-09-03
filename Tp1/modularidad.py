import networkx as nx
 
def modularidad(g, atributo):
    A = nx.adjacency_matrix(g)
    B = 0
    m = len(g.edges())
    k = []
    name = []
    for a, b in g.degree():
        k.append(b)
        name.append(a)
    for i in range(len(k)):
        for j in range(len(k)):
            if g.node[name[i]][atributo] == g.node[name[j]][atributo]:
                B += A[i, j] - k[i]*k[j]/(2*m)
    return float(B)/(2*m)

print modularidad(dolph, 'gender')
