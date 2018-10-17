# -*- coding: utf-8 -*-
"""
Created on Wed Oct 17 18:26:03 2018

@author: Gabo
"""

''' Definimos la funcion para leer el archivo geneX '''
import re
import numpy as np
import networkx as nx
import time
import sys
import matplotlib.pyplot as plt
import seaborn as sn
sys.path.append('./Tp3/')

from scipy.cluster.hierarchy import dendrogram
from scipy.cluster.hierarchy import linkage
from scipy.cluster.hierarchy import cophenet
from scipy.spatial.distance import pdist

def ldata(archive):
    f=open(archive)
    data=[]
    for line in f:
        line = re.sub('["\n]', '', line)
        col=line.split(',')
        data.append(col)
    return data
#%%
geneX=ldata('Tp3/geneX.csv')
valores_genes=np.zeros((500, 12))
nombre_genes=[]
for j in range(1,501):
    nombre_genes.append(geneX[j][0])
    for i in range(1,13):
        valores_genes[j-1,i-1]=float(geneX[j][i])
#%%
C = np.corrcoef(valores_genes) # Matriz de correlaciones
S = (1 + C) / 2 # Matriz de similaridades
A = np.array(S >= 0.95, dtype=int)

fig, (ax1, ax2) = plt.subplots(1, 2)
ax1.set_title('Matriz de similaridades')
sn.heatmap(S, ax=ax1, square=True)
ax2.set_title('Matriz de adyacencia')
sn.heatmap(A, ax=ax2, square=True, cbar=False)

#%%
fig, axes = plt.subplots(3,3, figsize=(9,9))
axes = np.ravel(axes)

umbrales = np.linspace(0, 1, 9)[::-1]
matrices = []
for ax, umbral in zip(axes, umbrales):
    arr = np.array(S >= umbral, dtype=int)
    matrices.append(arr)
    ax.set_title('Umbral = {}'.format(umbral))
    sn.heatmap(arr, ax=ax, cbar=False, square=True)
fig.tight_layout()

#%%
# Para la función linkage, podemos pasarle una matriz en la cual
# las filas tienen que ser vectores de observaciones, y las columnas
# tienen que corresponder a cada una de las coordenadas de estos vectores
# Cada gen es un vector en un espacio 12-dimensional

# Alternativamente, podemos pasarle un array 1d llamado "matriz de distancias
# condensada", que contiene los elementos de la submatriz triangular superior
# de la matriz de distancias

# Distancia == similaridad

def ij_to_idx(i, j, n):
    """
    Convierte índices (i,j) de cada entrada de una matriz n por n
    a posiciones en una "condensed distance matrix". Para esto,
    considera solo la submatriz triangular superior de la matriz
    original (i < j).
    
    This is very easy to understand:

    with i*n + j you go to the position in the square-formed matrix;
    with - i*(i+1)/2 you remove lower triangle (including diagonal) in all lines before i;
    with - i you remove positions in line i before the diagonal;
    with - 1 you remove positions in line i on the diagonal.

    Copiado de:
    https://stackoverflow.com/questions/13079563/how-does-condensed-distance-matrix-work-pdist#19581227
    """
    return i*n + j - i * (i + 1) // 2 - i - 1

def idx_to_ij_vectorized(idx, n):
    """
    Convierte array de índices idx a arrays de índices ii y jj de la matriz no
    condensada. Trabaja con matrices triangulares superiores al igual
    que la otra función.
    
    Copiado de:
    https://stackoverflow.com/questions/13079563/how-does-condensed-distance-matrix-work-pdist#19581227
    """
    n_row_elems = np.cumsum(np.arange(1, n)[::-1])
    ii = (n_row_elems[:, None] - 1 < idx[None, :]).sum(axis=0)
    shifts = np.concatenate([[0], n_row_elems])
    jj = np.arange(1, n)[ii] + idx - shifts[ii]
    return ii, jj
#%%
matriz_de_prueba = np.array([[1,2,3],
                             [4,5,6],
                             [7,8,9]])
n = 3
l = (n * (n - 1) // 2)    
cond_dist_prueba = np.zeros(l)
for i in range(n):
    for j in range(i + 1, n):
        cond_dist_prueba[ij_to_idx(i, j, n)] = matriz_de_prueba[i,j]
        
print(cond_dist_prueba)
# Éxito!

ii, jj = idx_to_ij_vectorized(np.arange(n * (n - 1) // 2), n)
print([matriz_de_prueba[i, j] for i, j in zip(ii, jj)])
# Bien, me da lo mismo y me va a resultar más útil
#%%
def matrix_to_condensed(matrix):
    n = matrix.shape[0]
    l = (n * (n - 1) // 2)
    assert matrix.shape == (n, n)
    condensed = np.zeros((l))
    idx = np.arange(l)
    ii, jj = idx_to_ij_vectorized(idx, n)
    for k in range(l):
        condensed[k] = matrix[ii[k], jj[k]]
    return condensed
#%%
    
# Convierto las similaridades a condensed
condensed_S = matrix_to_condensed(S)
linkage_matrix = linkage(condensed_S)
fig, ax = plt.subplots()
dendrogram(linkage_matrix)