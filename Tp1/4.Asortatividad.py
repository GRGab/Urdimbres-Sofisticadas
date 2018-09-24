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
degree_2, annd_2 = annd(net_science)
degree_1, annd_1 = annd(july)

#%%
#En este plot se puede ver que no aparece una recta, por lo que no hay que ajustar aca.

f, ([ax1, ax2], [ax3, ax4]) = plt.subplots(2, 2)
f.tight_layout()

plt.sca(ax1)
ax1.set_title('(a)  Lineal')
ax1.plot(degree_1, annd_1, '.')
ax1.set_ylabel(r'$k_{nn}$')
ax1.set_xlabel('k')

plt.sca(ax2)
ax2.set_title('(b)  Log Log')
ax2.loglog(degree_1, annd_1, '.')
ax2.set_ylabel(r'$k_{nn}$')
ax2.set_xlabel('k')
#ax2.yscale('log')
#ax2.xscale('log')

plt.sca(ax3)
coefs_reversed = np.flip(annd_1, 0)
cumulative = np.cumsum(coefs_reversed)
cumulative = np.flip(cumulative, 0)
ax3.set_title('(c)  Cumulative')
ax3.loglog(degree_1, cumulative, '.')
ax3.set_xlabel(r'$k_{nn}$')
#ax3.yscale('log')
#ax3.xscale('log')

plt.sca(ax4)
ax4.set_title('(d)  Log binned histogram')
histograma(annd_1, ax=ax4, xlabel=r'$k_{nn}$', labelsize=10, ticksize=10)


plt.savefig('Ej 4 plotteos del knn')
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
redes = [net_science, july, g_apms, g_y2h]
for G in redes:
    degree, annds = annd(G)
    log_k_nn = [np.log(i) for i in annds]
    log_k = [np.log(i) for i in degree]
    linear_model = Model(linear)
    data = RealData(log_k, log_k_nn)
    odr = ODR(data, linear_model, beta0=[0., 1.])
    out = odr.run()
    log_modelo = [j*out.beta[0]+out.beta[1] for j in log_k]
    
    plt.figure(redes.index(G))
    plt.plot(log_k, log_k_nn, '.')
    plt.plot(log_k, [i*out.beta[0]+out.beta[1] for i in log_k])
    #plt.plot(log_k[index_max], log_k_nn[index_max], 'o')
    #plt.plot(log_k[index_min], log_k_nn[index_min], 'o')
    r = nx.degree_assortativity_coefficient(G)
    
    print(r, out.beta[0])
    
'''
4.b
Si el grafico de k_nn vs k crece para k altos, significa que en promedio, los nodos con
grado alto van a tener vecinos que (en promedio) tambien tienen grado alto. Es decir,
los nodos de grado alto se van a conectar con nodos de grado alto.
Por otro lado, los nodos de k bajo tienen presentan k_nn bajo. Esto significa que los 
vecinos de estos nodos en promedio van a tener pocos vecinos. 
Luego, en este tipo de arreglos se observa que los nodos con grado alto se relacionan
con nodos de grado alto, y los nodos de grado bajo se relacionan con nodos de grado
bajo. En otras palabras, hubs se relacionan con hubs y los \textit{solitarios} con 
\textit{solitarios}. Es decir, la red es asortativa.
'''


