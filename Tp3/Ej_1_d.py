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
                xlabel=None, ylabel=None, titulo=None,
                labelsize_x=14, labelsize_y=18,
                sobrerrep=None):
    """
    Dibuja una matriz numérica como un heatmap de seaborn,
    asignándole un color a cada celda además de su número.

    El parámetro opcional sbrerrep es una matriz con enteros
    indicando de qué color se debe pintar cada celda. Sirve
    para matrices de p-valores del test exacto de Fisher,
    para indicar qué población es la que se ve sobrerrepresentada
    en cada cluster.
    """
    plt.figure(figsize=(8,6))
    df = pd.DataFrame(matriz, columns=columnas, index=filas)
    if sobrerrep is None:
        ax = sn.heatmap(df, fmt=".2f", vmin = 0, vmax=1, square=True,
                        annot=annot, annot_kws={'fontsize': 20})
    else:
        colores = []
        num_to_color = {-1: 'w', 0:'lightgrey', 1:'dodgerblue', 2:'r'}
        
        for num in np.ravel(sobrerrep, order='C'):
            colores.append(num_to_color[num])
        ax = sn.heatmap(df, fmt=".2f", square=True,
                        cbar=False, annot=annot,
                        annot_kws={'fontsize': 20, 'color': 'k'})
        ax.collections[0].set_facecolor(colores)
    ax.set_yticklabels(ax.get_yticklabels(), rotation = 0)
    ax.set_xticklabels(ax.get_xticklabels(), rotation = 0)
    for label in ax.get_yticklabels():
        label.set_size(labelsize_y)
    for label in ax.get_xticklabels():
        label.set_size(labelsize_x)
    if titulo is not None:
        ax.set_title(titulo, fontsize=max(labelsize_x, labelsize_y))
    if xlabel is not None:
        ax.set_xlabel(xlabel, fontsize=labelsize_x)
    if ylabel is not None:
        ax.set_ylabel(ylabel, fontsize=labelsize_y)
    plt.tight_layout()
    plt.show()

# Importamos particiones
dolph = read_gml('Tp3/dolphins.gml')    
lista_de_metodos = ["infomap","label_prop", "fastgreedy", "eigenvector",
                    "louvain", "edge_betweenness", "walktrap"]
nombres_metodos = ["Infomap","Label\nProp", "Fast\nGreedy", "Eigen", "Louvain",
                   "Edge\nBet", "Walktrap"]
# npzfile = np.load('Tp3/tc03Data/Ej_b_particiones_10000recableos.npz')
# mod_rewire = npzfile['mod_rewire']
# mod_original = npzfile['mod_original']
# sil_rewire = npzfile['sil_rewire']
# sil_original = npzfile['sil_original']

def calcular_fisher(particiones, genders, thres=0.05):
    """
    Para cada método y para cada cluster, evalúo la sobrerrepresentación
    de género a través del p-valor de un test exacto de Fisher
    OJO: tener en cuenta que no todos los nodos son macho o hembra

    Output
    ------
    `pfisher` : array 2d
        Matriz de p-valores de los tests
    `sobrerrep` : array 2d
        Cada entrada es 0 si p >= thres, y en caso contrario es 1 si
        hay más hembras que en la red total y 2 si hay menos.
    """
    n_metodos = len(lista_de_metodos)
    max_clusters = max(len(part) for part in particiones)
    n_hembras = sum(1 for key, val in genders.items() if val=='f')
    n_machos = sum(1 for key, val in genders.items() if val=='m')
    # En general vamos a ignorar a los delfines sin género asignado
    ntot_efectivo = n_hembras + n_machos
    # Inicializo la matriz de p-valores con y matriz de sobrrerepresentacxiones
    pfisher = np.full([max_clusters, n_metodos], np.nan)
    sobrerrep = np.full([max_clusters, n_metodos], np.nan)
    for i, particion in enumerate(particiones):
        for j, cluster in enumerate(particion):
            # La tabla tiene como filas ('en el cluster', 'no en el cluster')
            # y como columnas ('hembra', 'macho')
            # A los que no se sabe el género, simplemente los ignoramos
            # Es decir que efectivamente los sacamos del cluster antes de realizar
            # el análisis
            hembras_encluster = sum(1 for nodo in cluster if genders[nodo] == 'f')
            machos_encluster = sum(1 for nodo in cluster if genders[nodo] == 'm')
            table = [[hembras_encluster,             machos_encluster],
                    [n_hembras - hembras_encluster, n_machos - machos_encluster]]
            p = fisher_exact(table)[1]
            pfisher[j,i] = p
            # Rellenamos sobrerrep
            if pfisher[j,i] >= thres:
                sobrerrep[j,i] = 0
            else:
                nclust_efectivo = hembras_encluster + machos_encluster
                hembras_porazar = nclust_efectivo * n_hembras / ntot_efectivo
                if hembras_encluster > hembras_porazar:
                    sobrerrep[j, i] = 1
                else:
                    sobrerrep[j, i] = 2
    sobrerrep[np.isnan(sobrerrep)] = -1
    return pfisher, sobrerrep

particiones = np.load('Tp3/tc03Data/Ej_b_particiones_tomi.npz')['salida_grafo_original']
genders = dict(ldata('Tp1/tc01_data/dolphinsGender.txt'))
pfisher, sobrerrep = calcular_fisher(particiones, genders, thres=0.05)
plot_matrix(pfisher, columnas=nombres_metodos, titulo='p-valores (test de Fisher exacto)',
            ylabel='Comunidad')
plot_matrix(pfisher, columnas=nombres_metodos, titulo='p-valores (test de Fisher exacto)',
            ylabel='Comunidad', sobrerrep=sobrerrep)            
