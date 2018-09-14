import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from networkx.readwrite.gml import read_gml
from lectura import ldata
from histograma import histograma


net_science = read_gml('/home/matias/Documentos/Facultad/Redes/Urdimbres-Sofisticadas/Tp1/tc01_data/netscience.gml')
july = read_gml('/home/matias/Documentos/Facultad/Redes/Urdimbres-Sofisticadas/Tp1/tc01_data/as-22july06.gml')

 
lista_1 = nx.average_neighbor_degree(july)
lista_2 = nx.average_neighbor_degree(net_science)

vecinos_degree_1 = lista_1.values()
nodes_1 = lista_1.keys()
#    
vecinos_degree_2 = lista_2.values()
nodes_2 = lista_2.keys()
#%%
#k_1, values_1 = np.histogram((list(vecinos_degree_1)))
#k_2, values_2 = np.histogram(list(vecinos_degree_2))
#
#cumulative = np.cumsum(k)
## plot the cumulative function
#plt.figure()
#plt.plot(values_1[:-1], cumulative, '.b')
#plt.figure()
#plt.plot(values_2[:-1], cumulative, '.r')
#
#
##%%
#plt.figure()
#plt.plot(k, '.')
#plt.yscale('log')
#plt.xscale('log')
#
#plt.hist(vecinos_degree_1,bins=int(len(vecinos_degree_1)**0.5))
#
##%%
#sor
