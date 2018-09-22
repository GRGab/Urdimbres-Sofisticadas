from __future__ import division
import os
os.chdir('/home/tomas/Desktop/Redes complejas/Urdimbres-Sofisticadas/Tp1')
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
    """ Recibe una lista, y entrega el indice del elemento mas cercano a la unidad."""
    compliance = []
    for j in range(0, len(v)):
        compliance.append(abs(v[j] - 1))
    return compliance.index(np.min(compliance))

def nan_delete(k_nn):
    """Remueve los inf y los nan de la lista manteniendo su k especifico.
    .
    .
    """
    k = np.arange(0, len(k_nn), 1)
    k_nn_temp = []
    k_temp = []    
    for i in range(len(k_nn)):
        if not np.isinf(k_nn[i]) and not np.isnan(k_nn[i]) and not np.isinf(k[i]):
            k_nn_temp.append(k_nn[i])
            k_temp.append(k[i])
    return k_temp, k_nn_temp

def annd(red):
    """ annd (average neighbouhr degree distribution) devuelve el annd en orden
    de los grados del nodo.
    
    Returns: k, k_nn
    
    k: array con los grados de la red (eje X)
    k_nn: array con los annd promediados por grado (eje Y)
    .
    .
    """
    nombres = list(red.nodes)
    avnedeg = nx.average_neighbor_degree(red)
    grados = nx.degree(red)
    a = []
    for i in range(max(dict(nx.degree(red)).values())):
        b = []
        for j in range(len(nombres)):
            if i == grados[nombres[j]]:
                b.append(avnedeg[nombres[j]])
        a.append(np.mean(b))
    k, k_nn = nan_delete(a)
    return k, k_nn

def gamma(G, k, k_min):
    k_temp = [np.log(i/(k_min-0.5)) for i in k]
    return (1 + len(G.nodes) / np.sum(k_temp))


net_science = read_gml('/home/tomas/Desktop/Redes complejas/Urdimbres-Sofisticadas/Tp1/tc01_data/netscience.gml')
july = read_gml('/home/tomas/Desktop/Redes complejas/Urdimbres-Sofisticadas/Tp1/tc01_data/as-22july06.gml')

 
degree, annd = annd(july)


        
#En este plot se puede ver que no aparece una recta, por lo que no hay que ajustar aca.
plt.figure(1)
plt.plot(av_ne_degree_1, '.')
plt.figure(2)
plt.plot(av_ne_degree_2, '.')

#En este plot se puede ver que no aparece una recta, por lo que no hay que ajustar aca.
plt.figure(3)
plt.plot(degree, annd, '.')
plt.yscale('log')
plt.xscale('log')
plt.figure(4)
plt.plot(av_ne_degree_2, '.')
plt.yscale('log')
plt.xscale('log')


#En este plot se puede ver la cumulativa de los Knn en escala logaritmica.
coefs_reversed = np.flip(k_nn, 0)
cumulative = np.cumsum(coefs_reversed)
cumulative = np.flip(cumulative, 0)
plt.figure(5)
plt.plot(k, cumulative, '.')
plt.yscale('log')
plt.xscale('log')



#%%
#Metodo iterativo para calcular la mejor recta con el metodo del chi2

k_nn = [np.log(i) for i in annd]
k = [np.log(i) for i in degree]

def chi2_iterative(k, k_nn, Foward = True):
    """
    Esta funcion agarra un eje X (k), un eje Y (k_nn) y busca los parametros para
    ajustar la mejor recta, buscando el regimen lineal de la curva. Esto lo hace
    sacando puntos de la curva, ajustando la curva resultante, y luego comparando 
    los parametros de los distintos ajustes, seleccionando el de menor chi2.
    
    Si Foward=True entonces la funcion va a ir sacando puntos del final para
    encontrar kmax. Si Foward=False, la funcion va a sacar puntos del principio para
    calcular kmin. El punto va a estar dado por k[index].
    
    Returns: m, b, chi2_stat, index
    
    m: pendiente de la recta resultante
    b: ordenada de la recta resultante
    chi2: estadistico de chi2 de la recta resultante
    index: indice del elemento donde empieza/termina el regimen lineal.
    .
    .
    """
    chi2_list = []
    m_list = []
    b_list = []
    if Foward == True:
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
    else:
        for j in range(0, len(k)-3):
            k_nn_temp = k_nn[j:]
            k_temp = k[j:]
            linear_model = Model(linear)
            data = RealData(k_temp, k_nn_temp)
            odr = ODR(data, linear_model, beta0=[0., 1.])
            out = odr.run()
            chi2_list.append(out.res_var)
            m_list.append(out.beta[0])
            b_list.append(out.beta[1])        
    #index = ClosestToOne(chi2_list)
    index = chi2_list.index(min(chi2_list))
    m = m_list[index]
    b = b_list[index]
    chi2 = chi2_list[index]
    
    return m, b, chi2, index

m, b, ks_stat, index_max = chi2_iterative(k, k_nn, Foward=True)
k_max = np.exp(k[index_max])
_,_,_, index_min = chi2_iterative(k, k_nn, Foward=False)
k_min = np.exp(k[index_min])

plt.plot(k, k_nn, '.')
plt.plot(k, [i*m+b for i in k])
plt.plot(k[index_max], k_nn[index_max], 'o')
plt.plot(k[index_min], k_nn[index_min], 'o')

gamma = gamma(july, degree, k_min)


#%%
#Metodo iterativo para calcular la mejor recta con el metodo de kolmogorov smirnov

k_nn = [np.log(i) for i in annd]
k = [np.log(i) for i in degree]

def ks_iterative(k, k_nn, Foward = True):
    """
    Esta funcion agarra un eje X (k), un eje Y (k_nn) y busca los parametros para
    ajustar la mejor recta, buscando el regimen lineal de la curva. Esto lo hace
    sacando puntos de la curva, ajustando la curva resultante, y luego comparando 
    los parametros de los distintos ajustes con el metodo de Kolmogorov Smirnoff.
    
    Si Foward=True entonces la funcion va a ir sacando puntos del final para
    encontrar kmax. Si Foward=False, la funcion va a sacar puntos del principio para
    calcular kmin. El punto va a estar dado por k[index].
    
    Returns: m, b, ks_stat, index
    
    m: pendiente de la recta resultante
    b: ordenada de la recta resultante
    ks_stat: estadistico de KS de la recta resultante
    index: indice del elemento donde empieza/termina el regimen lineal.
    .
    .
    """
    KS_list = []
    pvalue_list = []
    m_list = []
    b_list = []
    
    if Foward==True:
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
    else:
        for j in range(0, len(k)-3):
            k_nn_temp = k_nn[j:]
            k_temp = k[j:]
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
    
    return m, b, ks_stat, index

m, b, ks_stat, index_max = ks_iterative(k, k_nn, Foward=True)
k_max = degree[index_max]
_,_,_, index_min = ks_iterative(k, k_nn, Foward=False)
k_min = degree[index_min]

plt.plot(k, k_nn, '.')
plt.plot(k, [i*m+b for i in k])
plt.plot(k[index_max], k_nn[index_max], 'o')
plt.plot(k[index_min], k_nn[index_min], 'o')

r = nx.degree_assortativity_coefficient(july)
gamma = gamma(july, degree, k_min)
