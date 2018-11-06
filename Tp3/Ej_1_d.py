import numpy as np
import matplotlib.pyplot as plt
plt.ion()
import pandas as pd
import seaborn as sn
from scipy.stats import fisher_exact
from lectura import ldata

import networkx as nx
from networkx import NetworkXError
from networkx.algorithms.community.community_utils import is_partition
from itertools import product
from networkx.readwrite.gml import read_gml
from histograma import histograma
import sys

def plot_matrix(matriz, filas=None, columnas=None, annot=True,
                labelsize_x=15, labelsize_y=20):
    df = pd.DataFrame(matriz, columns=columnas, index=filas)
    ax = sn.heatmap(df, fmt=".2f", vmin = 0, vmax=1, square=True,
                    annot=annot, annot_kws={'fontsize': 20} )
    ax.set_yticklabels(ax.get_yticklabels(), rotation = 0, fontsize = 12)
    ax.set_xticklabels(ax.get_xticklabels(), rotation = 0, fontsize = 12)
    for label in ax.get_yticklabels():
        label.set_size(labelsize_y)
    for label in ax.get_xticklabels():
        label.set_size(labelsize_x)
    plt.show()
plot_matrix(pfisher, columnas=nombres_metodos)

# Importamos particiones
dolph = read_gml('Tp3/dolphins.gml')    
lista_de_metodos = ["infomap","label_prop", "fastgreedy", "eigenvector",
                    "louvain", "edge_betweenness", "walktrap"]
nombres_metodos = ["Infomap","Label Prop", "Fastgreedy", "Eigenvectors", "Louvain",
                   "Edge Betweenness", "Walktrap"]
# npzfile = np.load('Tp3/tc03Data/Ej_b_particiones_10000recableos.npz')
# mod_rewire = npzfile['mod_rewire']
# mod_original = npzfile['mod_original']
# sil_rewire = npzfile['sil_rewire']
# sil_original = npzfile['sil_original']

particiones = np.load('Tp3/tc03Data/Ej_b_particiones_tomi.npz')['salida_grafo_original']
genders = dict(ldata('Tp1/tc01_data/dolphinsGender.txt'))

# Para cada método y para cada cluster, evalúo la sobrerrepresentación
# de género a través del p-valor de un test exacto de Fisher
n_metodos = len(lista_de_metodos)
max_clusters = max(len(part) for part in particiones)
n_hembras = sum(1 for key, val in genders.items() if val=='f')
n_machos = sum(1 for key, val in genders.items() if val=='m')
# OJO: tener en cuenta que no todos los nodos son macho o hembra

# Inicializo la matriz con nans
pfisher = np.full([max_clusters, n_metodos], np.nan)
for i, particion in enumerate(particiones):
    for j, cluster in enumerate(particion):
        # La tabla tiene como filas ('en el cluster', 'no en el cluster')
        # y como columnas ('hembra', 'macho')
        # A los que no se sabe el género, simplemente los ignoramos
        n_cluster = len(cluster)
        hembras_encluster = sum(1 for nodo in cluster if genders[nodo] == 'f')
        machos_encluster = sum(1 for nodo in cluster if genders[nodo] == 'm')
        table = [[hembras_encluster,             machos_encluster],
                 [n_hembras - hembras_encluster, n_machos - machos_encluster]]
        pfisher[j,i] = fisher_exact(table)[1]

plot_matrix(pfisher, columnas=nombres_metodos)