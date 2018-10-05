from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from lectura import ldata
import sys
sys.path.append('./Tp2/')
import pandas as pd
from scipy.odr import Model, RealData, ODR

ess = ldata('tc02Data/Essential_ORFs_paperHe.txt')
ap = ldata('tc02Data/yeast_AP-MS.txt')
lit = ldata('tc02Data/yeast_LIT.txt')
lit_r = ldata('tc02Data/yeast_LIT_Reguly.txt')
y2h = ldata('tc02Data/yeast_Y2H.txt')

g_ap = nx.Graph()
g_ap.add_edges_from(ap)
g_ap = max(nx.connected_component_subgraphs(g_ap), key=len)

g_lit = nx.Graph()
g_lit.add_edges_from(lit)
g_lit = max(nx.connected_component_subgraphs(g_lit), key=len)

g_y2h = nx.Graph()
g_y2h.add_edges_from(y2h)
g_y2h = max(nx.connected_component_subgraphs(g_y2h), key=len)

redes = ['AP', 'LIT', 'Y2H']
n_nodos = [len(g_ap.nodes()), len(g_lit.nodes()), len(g_y2h.nodes())]
n_enlaces = [len(g_ap.edges()), len(g_lit.edges()), len(g_y2h.edges())]
k_medio = [np.mean(dict(g_ap.degree).values()), np.mean(dict(g_lit.degree).values()), 
           np.mean(dict(g_y2h.degree).values())]
clustering_medio = [nx.average_clustering(g_ap), nx.average_clustering(g_lit), nx.average_clustering(g_y2h)]

d = {'Redes':redes, 'Nodos':n_nodos, 'Enlaces':n_enlaces, 'Grado medio':k_medio, 'Coef clustering medio': clustering_medio}
df = pd.DataFrame(data=d)
print(df)



#%% 
# Le agrego categoria de esencial a los nodos correspondientes.

def agregar_importancia(G):
    g = list(G.nodes())
    value = np.zeros([len(g)])
    for h in range(len(g)):
        for i in range(len(ess)):
            for j in range(len(ess[i])):
                if g[h] == ess[i][j]:
                    value[h] = 1
                    break
    return g, value



def ess_vs_k(g):
    nodos, values = agregar_importancia(g)
    nodos = []
    k = []
    for a, b in g.degree():
        nodos.append(a)
        k.append(b)
    k_norm = [i/max(k) for i in k]
    threshold = np.arange(0.01, 1, 0.001)
    y = []
    x = []
    for j in threshold:
        y_i = 0
        x_i = 0
        for i in range(len(nodos)):
            if k_norm[i]>j:
                if values[i] == 1:
                    y_i += 1
                x_i += 1
        if x_i != 0:
            y.append(y_i/x_i)
            x.append(x_i/len(nodos))
    return y, x

y_0, x_0 = ess_vs_k(g_ap)
y_1, x_1 = ess_vs_k(g_lit)
y_2, x_2 = ess_vs_k(g_y2h)

plt.plot(x_0, y_0, label = 'Ap')
plt.plot(x_1, y_1, label = 'Lit')
plt.plot(x_2, y_2, label = 'Y2H')
plt.legend()

#%%
# Ej D-1. Recontra chequear.

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

def desarme(g):
    nodos, values = agregar_importancia(g)
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


k, y_0 = desarme(g_lit)
y_1 = [np.log(1-i) for i in y_0]
plt.plot(k[:10], y_1[:10], 'o', label = 'Lit')
plt.legend()
linear_model = Model(Linear)
data = RealData(k[:10], y_1[:10])
odr = ODR(data, linear_model, beta0=[0., 1.])
out = odr.run()
plt.plot(k[:10], [out.beta[0]*k[i]+out.beta[1] for i in range(1,11)], 'b-')

k, y_0 = desarme(g_ap)
y_1 = [np.log(1-i) for i in y_0]
plt.plot(k[:10], y_1[:10], 'g^', label = 'Ap')
plt.legend()
linear_model = Model(Linear)
data = RealData(k[:10], y_1[:10])
odr = ODR(data, linear_model, beta0=[0., 1.])
out = odr.run()
plt.plot(k[:10], [out.beta[0]*k[i]+out.beta[1] for i in range(1,11)], 'g-')

k, y_0 = desarme(g_y2h)
y_1 = [np.log(1-i) for i in y_0]
plt.plot(k[:10], y_1[:10], 'r*', label = 'Y2H')
plt.legend()
linear_model = Model(Linear)
data = RealData(k[:10], y_1[:10])
odr = ODR(data, linear_model, beta0=[0., 1.])
out = odr.run()
plt.plot(k[:10], [out.beta[0]*k[i]+out.beta[1] for i in range(1,11)], 'r-')

plt.ylabel('log(1-P)')
plt.xlabel('Protein connectivity (k)')