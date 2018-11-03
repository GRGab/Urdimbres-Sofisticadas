#%%
import numpy as np
from networkx.readwrite.gml import read_gml
import numpy as np
import matplotlib.pyplot as plt
plt.ion()
import matplotlib.cm as cm
import networkx as nx
from networkx.algorithms.community.community_utils import is_partition
from networkx.readwrite.gml import read_gml
import sys
sys.path.append('./Tp3')
from funciones_tp3 import (calcular_particion, NotAPartition, indices_to_nodos_particion,
                           comunidad_a_color, histograma_recableo, histograma_recableo_multimetodo)
from guardar_particiones import guardar_particiones

from time import time

#%%
dolph = read_gml('Tp3/dolphins.gml')    
lista_de_metodos = ["infomap","label_prop", "fastgreedy", "eigenvector",
                    "louvain", "edge_betweenness", "walktrap"]
nombres_metodos = ["Infomap","Label Prop", "Fastgreedy", "Eigenvectors", "Louvain",
                   "Edge Betweenness", "Walktrap"]
#%%
#### Recablear y guardar
# t = time()
# guardar_particiones(dolph, 200, 10000, lista_de_metodos,
#                     output_path='Tp3/tc03Data/Ej_b_particiones_10000recableos',
#                     silencioso=True)
# print(time() - t, 'segundos')
# En compu escritorio Gabo tardó 50 minutos

#%%
#### Importar datos guardados
npzfile = np.load('Tp3/tc03Data/Ej_b_particiones_10000recableos.npz')
mod_rewire = npzfile['mod_rewire']
mod_original = npzfile['mod_original']
sil_rewire = npzfile['sil_rewire']
sil_original = npzfile['sil_original']

#%%
#### Análisis de modularidad
nbins = 15
# Histogramas por separado
for i, metodo in enumerate(nombres_metodos):
    fig, _ = histograma_recableo(mod_rewire[i], mod_original[i], metodo,
                                'Modularidad', bins=nbins)
    # uso 'lista_de_metodos[i]' en vez de 'metodo' para que el nombre
    # del archivo quede sin espacios en blanco.
    fig.savefig('Tp3/graficos/10000recableos/hist_mod_{}.png'.format(lista_de_metodos[i]))

# Histogramas todos juntos
fig, _ = histograma_recableo_multimetodo(mod_rewire, mod_original,
                                         nombres_metodos, 'Modularidad',
                                         bins=nbins)
fig.savefig('Tp3/graficos/10000recableos/hist_mod_todos.png')

#%%
#### Análisis de silhouettes
# 1) Curvas de silhouette para los grafos reales: 
# hecho en Ej_1_a.py.
# 2) Para los recableos: tiene sentido analizar silhouettes
# para todos los métodos? Veamos qué tan seguido se obtienen
# particiones con un único cluster.

for metodo, sils_metodo in zip(lista_de_metodos, sil_rewire):
    acc1 = 0
    acc2 = 0
    for sils_recableo in sils_metodo:
        acc1 += 1
        if len(sils_recableo) == 0:
            acc2 += 1
    print(acc2 / acc1, 'de las veces que se aplicó el método',
    metodo, 'se detectó una única comunidad.')

# Análisis 10000 recableos
# Vemos que para infomap se obtienen clusters únicos un 14% de
# las veces, para label_prop un 78% de las veces, y nunca para
# los demás métodos (o casi nunca: con eigenvectors ocurrió
# 9 veces de 10000).
# ----------------
# Decisión: hacemos análisis de silhouettes para infomap
# excluyendo los casos en los que hubo una única comunidad
# y reconociendo lo problemático de este hecho. No hacemos
# el análisis para label_prop porque no tiene sentido.

# Eliminamos entonces a labelprop de nuestras listas
nombres_metodos_s = [nombres_metodos[i] for i in range(len(nombres_metodos)) if i != 1]
lista_de_metodos_s = [lista_de_metodos[i] for i in range(len(lista_de_metodos)) if i != 1]
sil_original = [sil_original[i] for i in range(len(sil_original)) if i != 1]
sil_rewire = [sil_rewire[i] for i in range(len(sil_rewire)) if i != 1]

# 3) Análisis de silhouettes promedios
# Calculamos los valores de silhouette observados para el grafo real
sils_flat = [[s for cluster in sil_vals for s in cluster] for sil_vals in sil_original]
sils_prom_original = [sum(s)/len(s) for s in sils_flat]
# Calculamos los valores medios de silhouette para los recableos
sils_prom_rewire = []
for sils_por_metodo in sil_rewire:
    sprom_por_metodo = []
    for sils_por_recableo in sils_por_metodo:
        # Si me encuentro con un [] es porque no se calcularon
        # los coefs de silhouette en ese caso.
        # Por lo tanto chequeo que no sea un [] antes de promediar
        if len(sils_por_recableo) > 0:
            # Aplasto todos los valores a una única lista
            flattened = [s for cluster in sils_por_recableo for s in cluster]  
            # Promedio dicha lista
            sprom = sum(flattened) / len(flattened)
            sprom_por_metodo.append(sprom)
        else:
            # Si era un [], entonces agrego un nan
            sprom_por_metodo.append(np.nan)
    sils_prom_rewire.append(sprom_por_metodo)

# Elimino los nans antes de pasarle los valores
# a la función histograma.
for i in range(len(sils_prom_rewire)):
    arr = np.array(sils_prom_rewire[i])
    arr = arr[~np.isnan(arr)]
    sils_prom_rewire[i] = list(arr)

# Graficamos
nbins = 15
# Histogramas por separado
for i, metodo in enumerate(nombres_metodos_s):
    fig, _ = histograma_recableo(sils_prom_rewire[i],
                                 sils_prom_original[i],
                                 metodo, 'Silhouette medio', bins=nbins)
    # uso 'lista_de_metodos[i]' en vez de 'metodo' para que el nombre
    # del archivo quede sin espacios en blanco.
    fig.savefig('Tp3/graficos/10000recableos/hist_sprom_{}.png'.format(lista_de_metodos_s[i]))

# Histogramas todos juntos.
fig, _ = histograma_recableo_multimetodo(sils_prom_rewire, sils_prom_original,
                                         nombres_metodos_s, 'Silhouette medio',
                                         bins=nbins)
fig.savefig('Tp3/graficos/10000recableos/hist_sprom_todos.png')
