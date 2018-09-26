#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  9 21:45:45 2018

@author: tomas
"""

import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from networkx.readwrite.gml import read_gml
from histograma import histograma
from __future__ import division

internet = read_gml('Tp1/tc01_data/as-22july06.gml')
nodes = []
degrees = []
for a, b in internet.degree():
    nodes.append(a)
    degrees.append(b)
#%% Comparación de visualizaciones

# Para comparar los bineados logarítmicos y no logaasdrítmicos, lo justo es
# excluir a los nodos de grado 0 en ambos 

fig, axes = plt.subplots(4, 2, figsize=(10,10))
axes = axes.flatten()
logbinss = [0, 0, 0, 0, 1, 1, 1, 1]
logxs    = [0, 0, 1, 1, 0, 0, 1, 1]
logys    = [0, 1, 0, 1, 0, 1, 0, 1]

t = ['Bines lineales', 'Bines logarítmicos']
titulos  = [t[i] for i in logbinss]
xlabels = [('Grado (adim.)' if i in [6,7] else None) for i in range(8)]
ylabels = [(True if i % 2 == 0 else False) for i in range(8)]

for i in range(8):
    histograma(degrees,
               logbins=logbinss[i], logx=logxs[i], logy=logys[i], ax=axes[i],
               titulo=titulos[i], xlabel=xlabels[i], ylabel=ylabels[i],
               ecolor='k', errbars=False, 
               labelsize=10, ticksize=10,
               bins=(1, max(degrees) + 2, 100))
#%%
'''Algo interesante de ver es que el 98% de los degrees estan entre los 
degrees 0 y 20:'''
frac = np.sum([d <= 20 for d in degrees]) / len(degrees)
print(frac)
#%%

import rpy2.robjects as ro # Al hacer esto se inicializa un subproceso de R
from rpy2.robjects.packages import importr
# Usando importr, importamos paquetes de R que van a funcionar algo 
# así como módulos de Python

#%%
## EJECUTAR ESTO si no tienen instalado el paquete igraph (para instalarlo)
## import rpy2's package module
## select a mirror for R packages
#utils = importr('utils')
#utils.chooseCRANmirror(ind=2) # elijo de dónde descargar el paquete
## Instalo
#from rpy2.robjects.vectors import StrVector
#utils.install_packages(StrVector(['igraph']))
#%%
# Realizo el ajuste de la powerlaw
igraph = importr('igraph')
# Creamos un vector de R pasándole los degrees
degrees_r = ro.FloatVector(degrees)
# Documentación de fit_power_law:
# https://rdrr.io/cran/igraph/man/fit_power_law.html
resultado = igraph.fit_power_law(degrees_r, implementation='plfit')
print(resultado.r_repr())

#%%
# Graficamos histograma + ajuste
kmin = resultado.rx2('xmin')[0]
gamma = resultado.rx2('alpha')[0]
ksp = resultado.rx2('KS.p')[0]

from scipy.special import zeta
def powerlaw(x, gamma, kmin):
    # Como nuestro ajuste fue sobre una distribución de probabilidad discreta,
    # la cte de normalización es 1 sobre la función zeta de Riemann generalizada
    return x**(-gamma) / zeta(gamma, kmin)

# ACLARACION IMPORTANTE
# Para que se grafique bien la ley de potencias, es necesario llamar
# a la función powerlaw poniendo kmin=1. Esto se debe a que el histograma
# de grados está normalizado arrancando desde k=1, y no desde el kmin que
# elige la función fit_power_law.

fig, ax = plt.subplots(figsize=(8,6))
titulo = 'Histograma de grados con ajuste por ley de potencias'
histograma(degrees, logbins=True, ax=ax, titulo=titulo,
           logx=True, logy=True,
           xlabel='k (adim.)', ylabel=True, ecolor='k', errbars=False, 
           labelsize=18, ticksize=16, bins=(1, max(degrees) + 2, 50))
xs = np.linspace(1, max(degrees) + 2, 1000)
ax.plot(xs, powerlaw(xs, gamma, 1), '--', color='deeppink',
        label=r'$p(k) \propto k^{-\gamma}$')
#xs = np.arange(1, max(degrees) + 2)
#ax.plot(xs, powerlaw(xs, gamma, kmin), 'o', color='deeppink',
#        label=r'$\gamma = $' + '{:.4g}'.format(gamma))
ax.plot([], [], ' ', label=r'$\gamma = $' + '{:.4g}'.format(gamma))
ax.plot([], [], ' ', label=r'$K_{min} = $' + '{:.0f}'.format(kmin))
ax.plot([], [], ' ', label='p-value (KS) = {:.2g}'.format(ksp))
ax.legend()
#ax.axvline(kmin, color='deeppink')
