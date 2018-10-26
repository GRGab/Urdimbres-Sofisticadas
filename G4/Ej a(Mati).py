import networkx as nx
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import igraph as igraph
import time

from os.path import join as osjoin
import sys
sys.path.append('./Tp3/')
from histograma import histograma
from funciones_ej_a_mati import ldata, calculate_infomap, modularity

from scipy.cluster.hierarchy import dendrogram, linkage, cophenet
from scipy.spatial.distance import pdist

def genes_toarray(geneX):
    valores_genes=np.zeros((500, 12))
    nombre_genes=[]
    for j in range(1,501):
        nombre_genes.append(geneX[j][0])
        for i in range(1,13):
            valores_genes[j-1,i-1]=float(geneX[j][i])
    return valores_genes 

def crear_ma(s, umbral):
   ma =  np.array(s>=umbral, dtype = int)
   return ma

def atribuir_comunidades(criterio):
    lista = []
    for i in range (min(criterio),max(criterio)+1):
        lista_por_comunidades = np.where(info==i)[0].tolist()
        if len(lista_por_comunidades) > 0:
            lista.append(lista_por_comunidades)
    return lista

#%% Punto 1    
geneX = ldata('Tp3/geneX.csv')
valores_genes = genes_toarray(geneX)
#%%
#Ejercicio a
r  = np.corrcoef(valores_genes) 
s = (1+r)/2
ma = crear_ma(s, umbral=0.95)
#%%
fig, axes = plt.subplots(3,3, figsize=(9,9))
axes = np.ravel(axes)

#plt.figure(1)
#plt.imshow(r, interpolation='nearest', cmap=plt.cm.seismic_r, extent=(0.5,10.5,0.5,10.5))
#plt.colorbar()
#
#plt.figure(2)
#plt.imshow(s, interpolation='nearest', cmap=plt.cm.ocean, extent=(0.5,10.5,0.5,10.5))
#plt.colorbar()
#
#plt.show()
        
#%%Punto b
graph=nx.from_numpy_matrix(ma) 
#
info  = np.array(calculate_infomap(ma))
lista_infomap = atribuir_comunidades(info)
mod_infomap = modularity(graph,lista_infomap)

#fastgreed = np.array(calculate_infomap(ma, method='fastgreedy'))
#lista_fastgreed = atribuir_comunidades(fastgreed)
#mod_fastgreed = modularity(graph,lista_fastgreed)
#
#nx.draw(graph,  node_size = 10)
#%% Ejercicio 2
#Dist entre todos los elementos segun la metrica.
#matriz de dist(md): shortest path entre nodos (u otra), o la s
pdist = pdist(s, metric='euclidean')
#Como asocia las ramas para el dendograma deendiendo de la distancia entre nodos
lin=linkage(pdist)
d = dendrogram(lin, get_leaves=False)





