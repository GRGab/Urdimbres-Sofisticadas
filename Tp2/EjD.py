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
ess =  ess[2:-4]
ess = [fila[1] for fila in ess]
ess = np.unique(ess)

apms = ldata('Tp2/tc02Data/yeast_AP-MS.txt')
lit = ldata('Tp2/tc02Data/yeast_LIT.txt')
lit_r = ldata('Tp2/tc02Data/yeast_LIT_Reguly.txt')
lit_r = [fila[:2] for fila in lit_r[1:]]
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
    k_posibles = np.arange(max(k)+1)
    y = []
    x = []
    for j in k_posibles:
        y_i = 0
        x_i = 0 # nro de nodos con grado k
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
import statsmodels.api as sm
nombres = {g_apms: "APMS", g_lit: "Lit", g_lit_reg: "Lit_reg",
           g_y2h: "Y2H"}

with plt.style.context(('seaborn')):
    fig, ax = plt.subplots(2, 2, figsize=(12,8),
                           sharex=True)
ax = np.ravel(ax)

ms, bs = np.zeros((2, 4)) # 4 grafos, 2 parámetros a ajustar
sigma_ms, sigma_bs = np.zeros((2, 4))

grafos = [g_lit, g_lit_reg, g_apms, g_y2h]
for i, g in enumerate(grafos):
    k_0, y_0 = desarme(g, ess)
    k_1, y_1 = k_0[:10], y_0[:10]
    y_2 = [np.log(1 - yi) for yi in y_1]
    k_2 = sm.add_constant(k_1) # para que haya param ordenada al origen
    results = sm.OLS(y_2, k_2).fit()
    # El orden de los params es 1) ordenada al origen, 2) pendiente
    bs[i] = results.params[0]
    sigma_bs[i] = results.bse[0]
    ms[i] = results.params[1]
    sigma_ms[i] = results.bse[1]
    # print('Grafo: {}'.format(nombres[g]))
    # print(results.summary())
    
    # Dibujamos
    xs = np.array(k_1)
    ys = y_2
    fontsize = 18
    ticksize = 16
    ax[i].plot(xs, ys, 'o', color='dodgerblue')
    ax[i].plot(xs, ms[i] * xs + bs[i], '-', label='Ajuste lineal',
               color='deeppink')
    ax[i].tick_params(labelsize=ticksize)
    if i in [2, 3]:
        ax[i].set_xlabel('Grado', fontsize=fontsize)
    if i in [0, 2]:
        ax[i].set_ylabel(r'$\log(1-P_E)$', fontsize=fontsize)
    ax[i].set_title(nombres[g], fontsize=fontsize)
    
fig.tight_layout()
fig.savefig('Tp2/Ej d/figura2b_He.png')

# Obtenemos alfas y betas
alfas = 1 - np.exp(ms)
betas = 1 - np.exp(bs)
# Propagación de errores (chequeé con Monte Carlos que la propagación es correcta)
sigma_alfas = np.exp(ms) * sigma_ms
sigma_betas = np.exp(bs) * sigma_bs

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
alpha_lit_recta = 1 - np.exp(out.beta[0])
plt.plot(k[:10], [out.beta[0]*k[i]+out.beta[1] for i in range(1,11)], 'b-')
params_lit = out

k, y_0 = desarme(g_apms, ess)
y_1 = [np.log(1-i) for i in y_0]
plt.plot(k[:10], y_1[:10], 'g^', label = 'Ap')
plt.legend()
linear_model = Model(Linear)
data = RealData(k[:10], y_1[:10])
odr = ODR(data, linear_model, beta0=[0., 1.])
out = odr.run()
alpha_apms_recta = 1 - np.exp(out.beta[0])
plt.plot(k[:10], [out.beta[0]*k[i]+out.beta[1] for i in range(1,11)], 'g-')
params_apms = out

k, y_0 = desarme(g_y2h, ess)
y_1 = [np.log(1-i) for i in y_0]
plt.plot(k[:10], y_1[:10], 'r*', label = 'Y2H')
plt.legend()
linear_model = Model(Linear)
data = RealData(k[:10], y_1[:10])
odr = ODR(data, linear_model, beta0=[0., 1.])
out = odr.run()
alpha_y2h_recta = 1 - np.exp(out.beta[0]) 
plt.plot(k[:10], [out.beta[0]*k[i]+out.beta[1] for i in range(1,11)], 'r-')
params_y2h = out

k, y_0 = desarme(g_lit_reg, ess)
y_1 = [np.log(1-i) for i in y_0]
plt.plot(k[:10], y_1[:10], 'k.', label = 'Lit_reg')
plt.legend()
linear_model = Model(Linear)
data = RealData(k[:10], y_1[:10])
odr = ODR(data, linear_model, beta0=[0., 1.])
out = odr.run()
alpha_lit_reg_recta = 1 - np.exp(out.beta[0])
plt.plot(k[:10], [out.beta[0]*k[i]+out.beta[1] for i in range(1,11)], 'k-')
params_lit_reg = out

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

def vecinos_comunes(G, nodo1, nodo2):
    vecinos1 = set(G[nodo1])
    vecinos2 = set(G[nodo2])
    return len(vecinos1.intersection(vecinos2))

def calcular_pares(G, numvecinos):
    num_pares = 0
    num_pares_mismotipo = 0
    for n1, n2 in nx.non_edges(G):
        assert n1 not in G[n2]
        if vecinos_comunes(G, n1, n2) >= numvecinos:
            num_pares += 1
            if G.nodes()[n1]['esencialidad'] == G.nodes()[n2]['esencialidad']:
                num_pares_mismotipo += 1
    return num_pares, num_pares_mismotipo

#%%
# Para las 4 redes:
for nom, g, numvec in zip(['LIT', 'AP', 'Y2H', 'LIT_REG'],
                          [g_lit, g_apms, g_y2h, g_lit_reg],
                          [3, 3, 1, 3]):
    print('{}: '.format(nom), calcular_pares(g, numvec))
    
#AP:  (11569, 5875)
#LIT:  (718, 383)
#Y2H:  (23013, 15045)
#LIT_REG:  (10777, 6187)

# La red LIT_REG da casi lo mismo que en el paper (en el cual da
# 10777, 6143).

