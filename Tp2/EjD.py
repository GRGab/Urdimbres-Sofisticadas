from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from lectura import ldata

import sys
from time import time
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
rsquared_adj = np.zeros((4)) # R^2 corregidos

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
    rsquared_adj[i] = results.rsquared_adj
    print('Grafo: {}'.format(nombres[g]))
    print(results.summary())
    
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
#fig.savefig('Tp2/Ej d/figura2b_He.png')

# Obtenemos alfas y betas
alfas = 1 - np.exp(ms)
betas = 1 - np.exp(bs)
# Propagación de errores (chequeé con Monte Carlos que la propagación es correcta)
sigma_alfas = np.exp(ms) * sigma_ms
sigma_betas = np.exp(bs) * sigma_bs

#%%
# Generamos tabla de resultados

#ms_consigmas = ['{:.2f} +/- {:.2f}'.format(x, err) for x, err in zip(ms, sigma_ms)]
#bs_consigmas = ['{:.2f} +/- {:.2f}'.format(x, err) for x, err in zip(bs, sigma_bs)]
#alfas_consigmas = ['{:.2f} +/- {:.2f}'.format(x, err) for x, err in zip(alfas, sigma_alfas)]
#betas_consigmas = ['{:.2f} +/- {:.2f}'.format(x, err) for x, err in zip(betas, sigma_betas)]
#rsquared_prolijos = ['{:.2f}'.format(x) for x in rsquared_adj]
ms_consigmas = ['{:.3f} +/- {:.3f}'.format(x, err) for x, err in zip(ms, sigma_ms)]
bs_consigmas = ['{:.3f} +/- {:.3f}'.format(x, err) for x, err in zip(bs, sigma_bs)]
alfas_consigmas = ['{:.3f} +/- {:.3f}'.format(x, err) for x, err in zip(alfas, sigma_alfas)]
betas_consigmas = ['{:.3f} +/- {:.3f}'.format(x, err) for x, err in zip(betas, sigma_betas)]
rsquared_prolijos = ['{:.3f}'.format(x) for x in rsquared_adj]

tabla_ajustes = pd.DataFrame(data = {'Pendiente': ms_consigmas,
                                     'Ord. al or.': bs_consigmas,
                                     r'$R^2$'+'adj.': rsquared_prolijos,
                                     r'$\alpha$': alfas_consigmas,
                                     r'$\beta$': betas_consigmas},
                             index=['Lit', 'Lit_reg', 'APMS', 'Y2H'])
tabla_ajustes

#%%

# Implementación de Tomi (viejo)

# =============================================================================
# fig, ax = plt.subplots()
# 
# k, y_0 = desarme(g_lit, ess)
# y_1 = [np.log(1-i) for i in y_0]
# plt.plot(k[:10], y_1[:10], 'bo', label = 'Lit')
# plt.legend()
# linear_model = Model(Linear)
# data = RealData(k[:10], y_1[:10])
# odr = ODR(data, linear_model, beta0=[0., 1.])
# out = odr.run()
# alpha_lit_recta = 1 - np.exp(out.beta[0])
# plt.plot(k[:10], [out.beta[0]*k[i]+out.beta[1] for i in range(1,11)], 'b-')
# params_lit = out
# 
# k, y_0 = desarme(g_apms, ess)
# y_1 = [np.log(1-i) for i in y_0]
# plt.plot(k[:10], y_1[:10], 'g^', label = 'Ap')
# plt.legend()
# linear_model = Model(Linear)
# data = RealData(k[:10], y_1[:10])
# odr = ODR(data, linear_model, beta0=[0., 1.])
# out = odr.run()
# alpha_apms_recta = 1 - np.exp(out.beta[0])
# plt.plot(k[:10], [out.beta[0]*k[i]+out.beta[1] for i in range(1,11)], 'g-')
# params_apms = out
# 
# k, y_0 = desarme(g_y2h, ess)
# y_1 = [np.log(1-i) for i in y_0]
# plt.plot(k[:10], y_1[:10], 'r*', label = 'Y2H')
# plt.legend()
# linear_model = Model(Linear)
# data = RealData(k[:10], y_1[:10])
# odr = ODR(data, linear_model, beta0=[0., 1.])
# out = odr.run()
# alpha_y2h_recta = 1 - np.exp(out.beta[0]) 
# plt.plot(k[:10], [out.beta[0]*k[i]+out.beta[1] for i in range(1,11)], 'r-')
# params_y2h = out
# 
# k, y_0 = desarme(g_lit_reg, ess)
# y_1 = [np.log(1-i) for i in y_0]
# plt.plot(k[:10], y_1[:10], 'k.', label = 'Lit_reg')
# plt.legend()
# linear_model = Model(Linear)
# data = RealData(k[:10], y_1[:10])
# odr = ODR(data, linear_model, beta0=[0., 1.])
# out = odr.run()
# alpha_lit_reg_recta = 1 - np.exp(out.beta[0])
# plt.plot(k[:10], [out.beta[0]*k[i]+out.beta[1] for i in range(1,11)], 'k-')
# params_lit_reg = out
# 
# plt.ylabel('log(1-P)')
# plt.xlabel('Protein connectivity (k)')
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

def calcular_pares(G, numvecinos, guardar_grados=False):
    num_pares = 0
    num_pares_mismotipo = 0
    if guardar_grados:
        ks_pares = []
    for n1, n2 in nx.non_edges(G):
        if vecinos_comunes(G, n1, n2) >= numvecinos:
            num_pares += 1
            if G.nodes()[n1]['esencialidad'] == G.nodes()[n2]['esencialidad']:
                num_pares_mismotipo += 1
            if guardar_grados:
                k1, k2 = g.degree(n1), g.degree(n2)
                ks_pares.append([k1, k2])
    if guardar_grados:
        return num_pares, num_pares_mismotipo, ks_pares
    else:
        return num_pares, num_pares_mismotipo
#%%
grafos = [g_lit, g_lit_reg, g_apms, g_y2h] # Respetar este orden
pares, pares_mismotipo = np.zeros((2, 4))
grados_pares = []
numvecinos = [1, 3, 3, 1]

ti = time()
for i, g in enumerate(grafos):
    num_pares, num_pares_mismotipo, ks_pares = calcular_pares(g, numvecinos[i],
                                                              guardar_grados=True)
    pares[i] = num_pares
    pares_mismotipo[i] = num_pares_mismotipo
    grados_pares.append(np.array(ks_pares))
tf = time()
print(tf - ti, ' segundos')

def num_pares_esperados(g, ks, alfa, beta):
    """ks debe ser array de forma (n, 2) donde n es el número de pares de nodos
    no vecinos con 'numvecinos' o más vecinos comunes (o más en general, el
    conjunto de pares de nodos que satisface la condición deseada). Debe
    contener los grados de todos los nodos en esos pares.
    
    Devuelve
    --------
    rs : array
        probs de en un dado par, ambos nodos sean esenciales
    ss : array
        probs de en un dado par, ambos nodos sean no esenciales
    num_esperado : float
        el número de pares del mismo tipo esperado según el modelo de He"""
    k1 = ks[:,0]
    k2 = ks[:,1]
    rs = (1 - (1 - beta) * (1 - alfa)**k1) * (1 - (1 - beta) * (1 - alfa)**k2)
    ss = (1 - beta)**2 * (1 - alfa)**(k1 + k2)
    
    num_esperado = np.sum(rs) + np.sum(ss)
    return rs, ss, num_esperado

esperados = []
for i, g in enumerate(grafos):
    _, _, num_esperado = num_pares_esperados(g, grados_pares[i], alfas[i], betas[i])
    esperados.append(num_esperado)

#print('redes: g_lit, g_lit_reg, g_apms, g_y2h')
#print('pares: {:.0f}, {:.0f}, {:.0f}, {:.0f}'.format(*pares))
#print('pares_mismotipo: {:.0f}, {:.0f}, {:.0f}, {:.0f}'.format(*pares_mismotipo))
#print('esperados: {:.0f}, {:.0f}, {:.0f}, {:.0f}'.format(*esperados))

# Para tener incerteza en la estimación de los números esperados
from scipy.stats import norm

incertezas = []
zscores = []
n_historias = 1000
for i, g in enumerate(grafos):
    alfas_mc = norm.rvs(loc=alfas[i], scale=sigma_alfas[i], size=n_historias)
    betas_mc = norm.rvs(loc=betas[i], scale=sigma_betas[i], size=n_historias)
    resultados_mc = np.zeros((n_historias))
    for j in range(n_historias):
        _, _, num_esperado = num_pares_esperados(g, grados_pares[i],
                                                 alfas_mc[j], betas_mc[j])
        resultados_mc[j] = num_esperado
    incertezas.append(np.std(resultados_mc))
    zscores.append((pares_mismotipo[i] - esperados[i]) / incertezas[i])
    
print('redes: g_lit, g_lit_reg, g_apms, g_y2h')
print('pares: {:.0f}, {:.0f}, {:.0f}, {:.0f}'.format(*pares))
print('pares_mismotipo: {:.0f}, {:.0f}, {:.0f}, {:.0f}'.format(*pares_mismotipo))
print('esperados: {:.0f}, {:.0f}, {:.0f}, {:.0f}'.format(*esperados))
print('incertezas: {:.0f}, {:.0f}, {:.0f}, {:.0f}'.format(*incertezas))
print('Z-scores: {:.2f}, {:.2f}, {:.2f}, {:.2f}'.format(*zscores))

# Resultados con numvecinos = [3,3,3,3]:
#redes: g_lit, g_lit_reg, g_apms, g_y2h
#pares: 718, 10777, 11569, 522
#pares_mismotipo: 383, 6187, 5875, 352
#esperados: 378, 5791, 6074, 291
#incertezas: 18, 150, 3516, 60
#Z-scores: 0.27, 2.65, -0.06, 1.02


# La red LIT_REG da casi lo mismo que en el paper (en el cual da
# 10777, 6143).

# Resultados con numvecinos = [1,3,3,1]:
#redes: g_lit, g_lit_reg, g_apms, g_y2h
#pares: 9934, 10777, 11569, 23013
#pares_mismotipo: 5767, 6187, 5875, 15045
#esperados: 5008, 5791, 6074, 14204
#incertezas: 100, 144, 2505, 1737
#Z-scores: 7.56, 2.75, -0.08, 0.48

# Resultados con numvecinos = [1,1,1,1]:
#redes: g_lit, g_lit_reg, g_apms, g_y2h
#pares: 9934, 220167, 25915, 23013
#pares_mismotipo: 5767, 127589, 13004, 15045
#esperados: 5008, 118693, 12923, 14204
#incertezas: 83, 3007, 3884, 1822
#Z-scores: 9.16, 2.96, 0.02, 0.46

# En conclusión, se puede rechazar el modelo de He para las redes
# g_lit, g_lit_reg pero no para las redes apms e y2h.
#%%
# Hacemos la tabla prolija linda
pares_prolijo = [int(x) for x in pares]
pares_mismotipo_prolijo = [int(x) for x in pares_mismotipo]
esperados_prolijo = ['{} +/- {}'.format(int(x), int(y)) for x, y in zip(esperados, incertezas)]
zscores_prolijo = ['{:.2f}'.format(x) for x in zscores]

tabla5 = pd.DataFrame(data = {'Número total de pares': pares_prolijo,
                              '# pares del mismo tipo': pares_mismotipo_prolijo,
                              '# esperado de pares dle mismo tipo': esperados_prolijo,
                              'Z-scores': zscores_prolijo},
                      index=['Lit', 'Lit_reg', 'APMS', 'Y2H'])
tabla5

#%%
n = 10
p, pmt = np.zeros((2, n))
for i in range(n):
    p[i], pmt[i] = calcular_pares(g_apms, 3)
#%%
print('{:.2g} +/- {:.2g}'.format(np.average(p), np.std(p)))
print('{:.2g} +/- {:.2g}'.format(np.average(pmt), np.std(pmt)))