from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from lectura import ldata
import sys
sys.path.append('./Tp2/')
from agregar_esencialidad import agregar_esencialidad
import pandas as pd
from scipy.odr import Model, RealData, ODR

ess = ldata('Tp2/tc02Data/Essential_ORFs_paperHe.txt')
ess = [fila for fila in ess if (fila != [] and fila[0] != 'updated')]
ess = [fila[1] for fila in ess]

apms = ldata('Tp2/tc02Data/yeast_AP-MS.txt')
lit = ldata('Tp2/tc02Data/yeast_LIT.txt')
lit_r = ldata('Tp2/tc02Data/yeast_LIT_Reguly.txt')
lit_r = [fila[:2] for fila in lit_r]
y2h = ldata('Tp2/tc02Data/yeast_Y2H.txt')

g_apms = nx.Graph()
g_apms.add_edges_from(apms)
g_apms = max(nx.connected_component_subgraphs(g_apms), key=len)

g_lit = nx.Graph()
g_lit.add_edges_from(lit)
g_lit = max(nx.connected_component_subgraphs(g_lit), key=len)

g_y2h = nx.Graph()
g_y2h.add_edges_from(y2h)
g_y2h = max(nx.connected_component_subgraphs(g_y2h), key=len)

g_lit_reg = nx.Graph()
g_lit_reg.add_edges_from(lit_r)
g_lit_reg = max(nx.connected_component_subgraphs(g_lit_reg), key=len)

redes = ['AP', 'LIT', 'Y2H', 'LIT_REG']
n_nodos = [len(g_apms.nodes()),
           len(g_lit.nodes()),
           len(g_y2h.nodes()),
           len(g_lit_reg.nodes())]
n_enlaces = [len(g_apms.edges()),
             len(g_lit.edges()),
             len(g_y2h.edges()),
             len(g_lit_reg.edges())]
k_medio = [np.mean(list(dict(g_apms.degree).values())),
           np.mean(list(dict(g_lit.degree).values())), 
           np.mean(list(dict(g_y2h.degree).values())),
           np.mean(list(dict(g_lit_reg.degree).values()))]
clustering_medio = [nx.average_clustering(g_apms),
                    nx.average_clustering(g_lit),
                    nx.average_clustering(g_y2h),
                    nx.average_clustering(g_lit_reg)]

d = {'Redes':redes,
     'Nodos':n_nodos,
     'Enlaces':n_enlaces,
     'Grado medio':k_medio,
     'Coef clustering medio': clustering_medio}
df = pd.DataFrame(data=d)
print(df)


#%%
##### Ejercicio D: Esencialidad: Módulos biológicos vs. Interacciónes Esenciales

##### 1 Recontra chequear.

def Linear(M, x):
    """
    Funcion lineal para ajustar con el ODR:
    
    >>> from scipy.odr import Model, RealData, ODR
    >>> linear_model = Model(Linear)
    >>> data = RealData(X, Y, sx=X_err, sy=Y_err)
    >>> odr = ODR(data, linear_model, beta0=[0., 1.])
    >>> out = odr.run()
    
    >>> m = out.beta[0]
    >>> b = out.beta[1]
    >>> m_err = out.sd_beta[0]
    >>> b_err = out.sd_beta[1]        
    >>> chi2 = out.res_var
    .
    .
    """
    m, b = M
    return m*x + b

def desarme(g, ess):
    nodos, values = agregar_esencialidad(g, ess)
    nodos = []
    k = []
    for a, b in g.degree():
        nodos.append(a)
        k.append(b)
    threshold = np.arange(0, max(k), 1)
    y = []
    x = []
    for j in threshold:
        y_i = 0
        x_i = 0
        for i in range(len(nodos)):
            if k[i] == j:
                x_i += 1
                if values[i] == 1:
                    y_i += 1
        if x_i != 0:
            y.append(y_i/x_i)
            x.append(j)
    return x, y
#%%
fig, ax = plt.subplots()

k, y_0 = desarme(g_lit, ess)
y_1 = [np.log(1-i) for i in y_0]
plt.plot(k[:10], y_1[:10], 'bo', label = 'Lit')
plt.legend()
linear_model = Model(Linear)
data = RealData(k[:10], y_1[:10])
odr = ODR(data, linear_model, beta0=[0., 1.])
out = odr.run()
plt.plot(k[:10], [out.beta[0]*k[i]+out.beta[1] for i in range(1,11)], 'b-')

k, y_0 = desarme(g_apms, ess)
y_1 = [np.log(1-i) for i in y_0]
plt.plot(k[:10], y_1[:10], 'g^', label = 'Ap')
plt.legend()
linear_model = Model(Linear)
data = RealData(k[:10], y_1[:10])
odr = ODR(data, linear_model, beta0=[0., 1.])
out = odr.run()
plt.plot(k[:10], [out.beta[0]*k[i]+out.beta[1] for i in range(1,11)], 'g-')

k, y_0 = desarme(g_y2h, ess)
y_1 = [np.log(1-i) for i in y_0]
plt.plot(k[:10], y_1[:10], 'r*', label = 'Y2H')
plt.legend()
linear_model = Model(Linear)
data = RealData(k[:10], y_1[:10])
odr = ODR(data, linear_model, beta0=[0., 1.])
out = odr.run()
plt.plot(k[:10], [out.beta[0]*k[i]+out.beta[1] for i in range(1,11)], 'r-')

k, y_0 = desarme(g_lit_reg, ess)
y_1 = [np.log(1-i) for i in y_0]
plt.plot(k[:10], y_1[:10], 'k.', label = 'Lit_reg')
plt.legend()
linear_model = Model(Linear)
data = RealData(k[:10], y_1[:10])
odr = ODR(data, linear_model, beta0=[0., 1.])
out = odr.run()
plt.plot(k[:10], [out.beta[0]*k[i]+out.beta[1] for i in range(1,11)], 'k-')

plt.ylabel('log(1-P)')
plt.xlabel('Protein connectivity (k)')

# =============================================================================
# # El fit de g_lit_reg está dando horrible, revisar
# =============================================================================

# Faltaría implementar un chi cuadrado para ver bondad de ajuste

#%%
##### 2: Tabla 5 de Zotenko
##### Difference between the observed and expected number of pairs where
##### both proteins are either essential or nonessential.

"""Epígrafe de la tabla:

The total number of pairs refers to the number of nonadjacent protein pairs
with three or more common neighbors in the network. (Due to the sparsity of
the Y2H network, the statistics are calculated for nonadjacent pairs having
one or more neighbors in common.) The nodes in the pair are of ‘‘the same
type’’ if they are both essential or both nonessential."""

from agregar_esencialidad import agregar_esencialidad_dict

for g in [g_lit, g_apms, g_y2h, g_lit_reg]:
    agregar_esencialidad_dict(g, ess)
#%%
def calcular_pares_mala(G, numvecinos):
    acc = 0
    n = G.order()
    A = nx.adj_matrix(G)
    for i in range(n):
        for j in range(n_nodos):
            vecinos_comunes = A[i].multiply(A[j]) # and lógico element-wise
            num_vec_com = vecinos_comunes.sum()
#            import pdb; pdb.set_trace(); 
            if A[i,j] != 0 and num_vec_com >= numvecinos:
                acc += 1
                if 
    return acc/2

# Tarda muchísimo
# calcular_pares(g_lit, 3) -> 1085.5
# Está funcionando mal porque debería dar un número entero
#%%

def vecinos_comunes(G, nodo1, nodo2):
    vecinos1 = set(G[nodo1])
    vecinos2 = set(G[nodo2])
    return len(vecinos1.intersection(vecinos2))

def calcular_pares(G, numvecinos):
    acc1 = 0
    acc2 = 0
    n = G.order()
    for n1, d1 in dict(G.nodes).items(): # nodo 1, dict 1
        for n2, d2 in dict(G.nodes).items(): # nodo 2, dict 2
            if (n2 not in G[n1] and
                vecinos_comunes(G, n1, n2) >= numvecinos):
                acc1 += 1
#                import pdb; pdb.set_trace(); 
                if d1['esencialidad'] == d2['esencialidad']:
                    acc2 += 1
    num_pares = acc1 / 2
    num_pares_mismotipo = acc2 / 2
    return num_pares, num_pares_mismotipo

# Tarda menos y sí da entero!
# calcular_pares(g_lit, 3) -> (1020.0, 685.0)

# Para las 4 redes:
for nom, g, numvec in zip(['AP', 'LIT', 'Y2H', 'LIT_REG'],
                          [g_lit, g_apms, g_y2h, g_lit_reg],
                          [3, 3, 1, 3]):
    print('{}: '.format(nom), calcular_pares(g, numvec))
    
# AP:  (1020.0, 685.0)
# LIT:  (11982.5, 6288.5)
# Y2H:  (23754.5, 15786.5)
# LIT_REG:  (11514.0, 6924.0)

# Posibles problemas
#-------------------
    
# 1) LIT y Y2H no dan un número entero, hay que revisasr por qué.
# 2) Los números no son nada que ver con los del paper! Igual me parece que
# ninguna de nuestras redes es igual a ninguna de las del paper, con lo cual
# esto no sería un problema en sí mismo
    