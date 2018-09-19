import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from networkx.readwrite.gml import read_gml
from lectura import ldata
from histograma import histograma
from scipy.odr import Model, RealData, ODR
from scipy import stats

def linear(M, x):
    m, b = M
    return x*m + b

def ClosestToOne(v):
    compliance = []
    for j in range(0, len(v)):
        compliance.append(abs(v[j] - 1))
    return compliance.index(np.min(compliance))

def inf_delete(k, k_nn):
    '''Remueve los inf y los nan de la lista manteniendo su k especifico'''
    k_nn_temp = []
    k_temp = []    
    for i in range(len(k_nn)):
        if not np.isinf(k_nn[i]) and not np.isnan(k_nn[i]) and not np.isinf(k[i]):
            k_nn_temp.append(k_nn[i])
            k_temp.append(k[i])
    return k_temp, k_nn_temp


net_science = read_gml('/home/tomas/Desktop/Redes complejas/Urdimbres-Sofisticadas/Tp1/tc01_data/netscience.gml')
july = read_gml('/home/tomas/Desktop/Redes complejas/Urdimbres-Sofisticadas/Tp1/tc01_data/as-22july06.gml')

 
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
#Metodo iterativo para calcular la mejor recta con el metodo del chi2

coefs = av_ne_degree_2

#k = np.delete(coefs, 0)
#k_nn = np.log(coefs)
#k = np.log(np.arange(1, len(coefs)+1, 1))

#k_nn = np.log(coefs[np.where(coefs)[0]])
#k = np.log(np.arange(1, len(k_nn)+1))

k = np.arange(0, len(coefs), 1)
k_nn = [np.log(i) for i in coefs]
k = [np.log(i) for i in k]

k, k_nn = inf_delete(k, k_nn)

chi2_list = []
m_list = []
b_list = []
for j in range(0, len(k)-3):
    k_nn_temp = k_nn[:len(k_nn)-j]
    k_temp = k[:len(k)-j]
    linear_model = Model(linear)
    data = RealData(k_temp, k_nn_temp)
    odr = ODR(data, linear_model, beta0=[0., 1.])
    out = odr.run()
    chi2_list.append(out.res_var)
    m_list.append(out.beta[0])
    b_list.append(out.beta[1])

index = ClosestToOne(chi2_list)
m = m_list[index]
b = b_list[index]
chi2 = chi2_list[index]

plt.plot(k, k_nn, '.')
plt.plot(k, [i*m+b for i in k])

#%%
#Metodo iterativo para calcular la mejor recta con el metodo de kolmogorov smirnov

KS_list = []
pvalue_list = []
m_list = []
b_list = []
for j in range(0, len(k)-3):
    k_nn_temp = k_nn[:len(k_nn)-j]
    k_temp = k[:len(k)-j]
    linear_model = Model(linear)
    data = RealData(k_temp, k_nn_temp)
    odr = ODR(data, linear_model, beta0=[0., 1.])
    out = odr.run()
    modelo = [j*out.beta[0]+out.beta[1] for j in k_temp]
    KS_list.append(stats.ks_2samp(k_nn_temp, modelo)[0])
    pvalue_list.append(stats.ks_2samp(k_nn_temp, modelo)[1])
    m_list.append(out.beta[0])
    b_list.append(out.beta[1])

index = KS_list.index(min(KS_list))
m = m_list[index]
b = b_list[index]
ks_stat = KS_list[index]

plt.plot(k, k_nn, '.')
plt.plot(k, [i*m+b for i in k])