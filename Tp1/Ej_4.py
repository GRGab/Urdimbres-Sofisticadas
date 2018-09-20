import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from networkx.readwrite.gml import read_gml
from lectura import ldata
from histograma import histograma
from scipy.odr import Model, RealData, ODR

net_science = read_gml('/home/matias/Documentos/Facultad/Redes/Urdimbres-Sofisticadas/Tp1/tc01_data/netscience.gml')
july = read_gml('/home/matias/Documentos/Facultad/Redes/Urdimbres-Sofisticadas/Tp1/tc01_data/as-22july06.gml')

 
lista_1 = nx.average_neighbor_degree(july)
lista_2 = nx.average_neighbor_degree(net_science)

vecinos_degree_1 = list(lista_1.values())
nodes_1 = list(lista_1.keys())
    
vecinos_degree_2 = list(lista_2.values())
nodes_2 = list(lista_2.keys())
    

av_ne_degree_1, bines_1 = np.histogram(vecinos_degree_1, bins = np.arange(0, len(vecinos_degree_1), 1))
av_ne_degree_2, bines_2 = np.histogram(vecinos_degree_2, bins = np.arange(0, len(vecinos_degree_2), 1))




def annd(red):
    ''' annd (average neighbouhr degree distribution) devuelve el and en orden
    de los grados del nodo'''
    nombres = list(red.nodes)
    avnedeg = nx.average_neighbor_degree(red)
    grados = nx.degree(red)
#    a = {}
#    for nombre in nombres:
#        key = grados[nombre]
#        if key not in a.keys():
#            a[keys] = []
#        a[keys].append(avnedeg[nombre])
    a = []
    for grado in range (max(dict(nx.degree(red)).values())):
        for nodos in 
        b = []
        b.append(avnedeg[grado])
        a.append(avnedeg[grado])
        
        


for i in range ()


#En este plot se puede ver que no aparece una recta, por lo que no hay que ajustar aca.
plt.figure(1)
plt.plot(av_ne_degree_1, '.')
plt.figure(2)
plt.plot(av_ne_degree_2, '.')

#En este plot se puede ver que no aparece una recta, por lo que no hay que ajustar aca.
plt.figure(3)
plt.plot(av_ne_degree_1, '.')
plt.yscale('log')
plt.xscale('log')
plt.figure(4)
plt.plot(av_ne_degree_2, '.')
plt.yscale('log')
plt.xscale('log')


#En este plot se puede ver la cumulativa de los Knn en escala logaritmica.
coefs_reversed = np.flip(av_ne_degree_1, 0)
cumulative = np.cumsum(coefs_reversed)
cumulative = np.flip(cumulative, 0)
plt.figure(5)
plt.plot(cumulative, '.')
plt.yscale('log')
plt.xscale('log')



#%%
#Metodo iterativo para calcular la mejor recta.

coefs = av_ne_degree_1
def linear(x, m, b):
    return x*m + b

def ClosestToOne(v):
    compliance = []
    for j in range(0, len(v)):
        compliance.append(abs(v[j] - 1))
    return compliance.index(np.min(compliance))



#k = np.delete(coefs, 0)
#k_nn = np.log(coefs)
#k = np.log(np.arange(1, len(coefs)+1, 1))

k_nn = np.log(coefs[np.where(coefs)[0]])
k = np.log(np.arange(1, len(k_nn)+1))

chi2 = []
Rq_err_temp = []
m = []
b = []
for j in range(0, len(k-10)):
    
    
    
    
    np.delete(k_nn, len(k_nn) - j-1)
    k = np.delete(k, len(k) - j-1)
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
