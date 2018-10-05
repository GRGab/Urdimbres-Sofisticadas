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
apms = ldata('Tp2/tc02Data/yeast_AP-MS.txt')
lit = ldata('Tp2/tc02Data/yeast_LIT.txt')
lit_r = ldata('Tp2/tc02Data/yeast_LIT_Reguly.txt')
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

redes = ['AP', 'LIT', 'Y2H']
n_nodos = [len(g_apms.nodes()), len(g_lit.nodes()), len(g_y2h.nodes())]
n_enlaces = [len(g_apms.edges()), len(g_lit.edges()), len(g_y2h.edges())]
k_medio = [np.mean(list(dict(g_apms.degree).values())),
           np.mean(list(dict(g_lit.degree).values())), 
           np.mean(list(dict(g_y2h.degree).values()))]
clustering_medio = [nx.average_clustering(g_apms), nx.average_clustering(g_lit), nx.average_clustering(g_y2h)]

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

k, y_0 = desarme(g_lit, ess)
y_1 = [np.log(1-i) for i in y_0]
plt.plot(k[:10], y_1[:10], 'o', label = 'Lit')
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

plt.ylabel('log(1-P)')
plt.xlabel('Protein connectivity (k)')