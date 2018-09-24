from __future__ import division

import sys
sys.path.append('./Tp1/')
#import os
#os.chdir('/home/tomas/Desktop/Redes complejas/Urdimbres-Sofisticadas/Tp1')
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from networkx.readwrite.gml import read_gml
from lectura import ldata
from histograma import histograma
from scipy.odr import Model, RealData, ODR
from scipy import stats
from asortatividad_funciones import (linear, ClosestToOne, annd, gamma,
                                     chi2_iterative, ks_iterative)

#%%Importamos las redes
#Mati: Tomi te comente esto porque solo funciona en tu compu, me parece que los
#path que puse aca abajo funcionan en las compus de todes.

#net_science = read_gml('/home/tomas/Desktop/Redes complejas/Urdimbres-Sofisticadas/Tp1/tc01_data/netscience.gml')
#july = read_gml('/home/tomas/Desktop/Redes complejas/Urdimbres-Sofisticadas/Tp1/tc01_data/as-22july06.gml')

net_science = read_gml('tc01_data/netscience.gml')
july = read_gml('tc01_data/as-22july06.gml')
degree_1, annd_1 = annd(net_science)
degree_2, annd_2 = annd(july)
degree_3, annd_3 = annd(g_apms)
degree_4, annd_4 = annd(g_y2h)
#%%
#En este plot se puede ver que no aparece una recta, por lo que no hay que ajustar aca.
plt.figure(1)
plt.plot(degree_1, annd_1, '.')
plt.figure(2)
plt.plot(degree_2, annd_2, '.')

#En este plot se puede ver que no aparece una recta, por lo que no hay que ajustar aca.
plt.figure(3)
plt.plot(degree_1, annd_1, '.')
plt.yscale('log')
plt.xscale('log')

plt.figure(4)
plt.plot(degree_2, annd_2, '.')
plt.yscale('log')
plt.xscale('log')

#%%
#Para July
log_k_nn = [np.log(i) for i in annd_4]
log_k = [np.log(i) for i in degree_4]
#%%
#En este plot se puede ver la cumulativa de los Knn en escala logaritmica.
coefs_reversed = np.flip(k_nn, 0)
cumulative = np.cumsum(coefs_reversed)
cumulative = np.flip(cumulative, 0)
plt.figure(5)
plt.plot(k, cumulative, '.')
plt.yscale('log')
plt.xscale('log')
#%%
#Calculo del ajuste por chi2 y ploteo del mismo
m, b, ks_stat, index_max = chi2_iterative(log_k, log_k_nn, Foward=True)
k_max_chi = np.exp(k[index_max])
_,_,_, index_min = chi2_iterative(log_k, log_k_nn, Foward=False)
k_min_chi = np.exp(k[index_min])

plt.plot(k, k_nn, '.')
plt.plot(k, [i*m+b for i in k])
plt.plot(k[index_max], k_nn[index_max], 'o')
plt.plot(k[index_min], k_nn[index_min], 'o')

#%% 4.iii
# Hacer esto para las 4 redes!

#Calculo del ajuste por Kolmogorov y ploteo del mismo
#_, _, _, index_max = ks_iterative(log_k, log_k_nn, Foward=True)
#k_max_kol = degree_2[index_max]
#m,b,ks_stat, index_min = ks_iterative(log_k[:index_max], log_k_nn[:index_max], Foward=False)
##m,b,ks_stat, index_min = ks_iterative(log_k, log_k_nn, Foward=True)
#k_min_kol = degree_2[index_min]

linear_model = Model(linear)
data = RealData(log_k, log_k_nn)
odr = ODR(data, linear_model, beta0=[0., 1.])
out = odr.run()
log_modelo = [j*out.beta[0]+out.beta[1] for j in log_k]

plt.plot(log_k, log_k_nn, '.')
plt.plot(log_k, [i*out.beta[0]+out.beta[1] for i in log_k])
#plt.plot(log_k[index_max], log_k_nn[index_max], 'o')
#plt.plot(log_k[index_min], log_k_nn[index_min], 'o')
r = nx.degree_assortativity_coefficient(g_y2h)

print(r, out.beta[0])
#%%

gamma_kol = gamma(july, degree_2, k_min_kol)
gamma_chi = gamma(july, degree_2, k_min_chi)

