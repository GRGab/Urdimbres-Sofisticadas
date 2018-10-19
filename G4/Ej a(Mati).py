from lectura import ldata
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
from funciones_ej_a import calculate_infomap


genex = ldata('Tp3/geneX.csv')


valores_genes=np.zeros((500, 12))
nombre_genes=[]
for j in range(1,501):
    nombre_genes.append(geneX[j][0])
    for i in range(1,13):
        valores_genes[j-1,i-1]=float(geneX[j][i])

#%%
r  = np.corrcoef(valores_genes) 
s = (1+r)/2
#plt.figure(1)
#plt.imshow(r, interpolation='nearest', cmap=plt.cm.seismic_r, extent=(0.5,10.5,0.5,10.5))
#plt.colorbar()
#
#plt.figure(2)
#plt.imshow(s, interpolation='nearest', cmap=plt.cm.ocean, extent=(0.5,10.5,0.5,10.5))
#plt.colorbar()
#
#plt.show()
#        

def crear_ma(s, umbral=0.95):
   ma =  np.array(s<umbral, dtype = int)
   return ma

ma = crear_ma(s, umbral=0.5)
#for i in (0,0.95):
#    fig, ax = plt.subplots()
    
info  = np.array(calculate_infomap(ma))
