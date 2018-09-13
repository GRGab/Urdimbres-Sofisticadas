import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from networkx.readwrite.gml import read_gml
from lectura import ldata
from histograma import histograma
from scipy.odr import Model, RealData, ODR

net_science = read_gml('/home/tomas/Desktop/Redes complejas/Urdimbres-Sofisticadas/Tp1/tc01_data/netscience.gml')
july = read_gml('/home/tomas/Desktop/Redes complejas/Urdimbres-Sofisticadas/Tp1/tc01_data/as-22july06.gml')

 
lista = nx.average_neighbor_degree(net_science)

vecinos_degree = lista.values()
nodes = lista.keys()
    

coefs, values = np.histogram(vecinos_degree, bins = np.arange(0, len(vecinos_degree), 1))


#En este plot se puede ver que no aparece una recta, por lo que no hay que ajustar aca.
plt.figure(1)
plt.plot(coefs, '.')
plt.yscale('log')
plt.xscale('log')


#En este plot se puede ver la cumulativa de los Knn en escala logaritmica.
coefs_reversed = np.flip(coefs, 0)
cumulative = np.cumsum(coefs_reversed)
cumulative = np.flip(cumulative, 0)
plt.figure(2)
plt.plot(cumulative, '.')
plt.yscale('log')
plt.xscale('log')



#%%
#Metodo iterativo para calcular la mejor recta.

def linear(x, m, b):
    return x*m + b

def ClosestToOne(v):
    compliance = []
    for j in range(0, len(v)):
        compliance.append(abs(v[j] - 1))
    return compliance.index(np.min(compliance))

k = np.delete(coefs, 0)
k_nn = np.log(coefs)
k = np.log(np.arange(1, len(coefs)+1, 1))

chi2 = []
Rq_err_temp = []
m = []
b = []
for j in range(0, len(k-10)):
    k_nn = np.delete(k_nn, len(k_nn) - j)
    k = np.delete(k_nn, len(k) - j)
    linear_model = Model(linear)
    data = RealData(k, k_nn)
    odr = ODR(data, linear_model, beta0=[0., 1.])
    out = odr.run()
    chi2.append(out.res_var)
    m.append(out.beta[0])
    b.append(out.beta[1])

index = ClosestToOne(chi2)
m = m[index]
b = b[index]
chi2 = chi2[index]
