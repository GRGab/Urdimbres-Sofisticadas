import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from networkx.readwrite.gml import read_gml
from lectura import ldata
from histograma import histograma


net_science = read_gml('/home/tomas/Desktop/Redes complejas/Urdimbres-Sofisticadas/Tp1/tc01_data/netscience.gml')
july = read_gml('/home/tomas/Desktop/Redes complejas/Urdimbres-Sofisticadas/Tp1/tc01_data/as-22july06.gml')

 
lista = nx.average_neighbor_degree(july)

vecinos_degree = lista.values()
nodes = lista.keys()
    

k, values = np.histogram(vecinos_degree)

plt.figure(1)
plt.plot(k, '.')
plt.yscale('log')
plt.xscale('log')

plt.hist(vecinos_degree,bins=int(len(vecinos_degree)**0.5))
