import numpy as np
import networkx as nx
from networkx.readwrite.gml import read_gml
from lectura import ldata

import sys
sys.path.append('./Tp3/')
from funciones_tp3 import calcular_particion, comunidad_a_color
from silhouettes import silhouettes, graficar_silhouettes
import matplotlib.pyplot as plt
plt.ion()

#%%

######## EJ 1 A ... Visualizacion de las distintas particiones.
######## Gráficos coloreados según particiones

## Dado que las curvas de silhouette también requieren
## un coloreo, las creamos al mismo tiempo para facilitar
## la consistencia entre el mismo y el coloreo de los nodos.
## Además, queremos dibujar el grafo coloreado según

# Importamos el grafo y las particiones obtenidas según los distintos
# métodos. Agregamos información sobre género
dolph = read_gml('Tp3/dolphins.gml')
genders = dict(ldata('Tp1/tc01_data/dolphinsGender.txt'))
metodos = ["Infomap","Label Prop", "Fastgreedy", "Eigenvectors", "Louvain",
           "Edge Betweenness", "Walktrap"]
for nodo, dict_nodo in dict(dolph.nodes).items():
    dict_nodo['gender'] = genders[nodo] # agrego el sexo del delfín a su dict

particiones = np.load('Tp3/tc03Data/Ej_b_particiones_tomi.npz')['salida_grafo_original']

#### Parámetros comunes a todos los dibujos
# Tamaño de los nodos
ns = 35
# Layout
pos = nx.spring_layout(dolph)

#### Colores correspondientes a cada cluster y a los géneros
colores_nodos, colores_clusters = [], []
for particion in particiones:
    a, b = comunidad_a_color(dolph, particion)
    colores_nodos.append(a)
    colores_clusters.append(b)

def genero_a_color(gender):
    if gender=='m':
        return 'red'
    elif gender=='f':
        return 'dodgerblue'
    else:
        return 'green'
colores_por_genero = [genero_a_color(g) for g in nx.get_node_attributes(dolph, "gender").values()]    

#### Gráfico coloreado por género
fig, ax = plt.subplots()
nx.draw(dolph, node_color=colores_por_genero, pos=pos, node_size=ns)
plt.tight_layout()
fig.savefig('Tp3/graficos/silhouette_corregido/grafo_generos.png')

#### Gráficos coloreados por clusterización
fig, axes = plt.subplots(2, 4, figsize=(12,6))
axes = np.ravel(axes)
for i, metodo in enumerate(metodos):    
    plt.sca(axes[i])
    axes[i].set_title(metodo)
    nx.draw(dolph, node_size=ns, node_color=colores_nodos[i], pos=pos)
    # También guardamos el gráfico por separado
    fig_temp, ax_temp = plt.subplots()
    ax_temp.set_title(metodo)
    nx.draw(dolph, node_size=ns, node_color=colores_nodos[i], pos=pos)
    fig_temp.tight_layout()
    fig_temp.savefig('Tp3/graficos/silhouette_corregido/Grafo {}.png'.format(metodo))
axes[7].axis('off')
fig.tight_layout()
fig.savefig('Tp3/graficos/silhouette_corregido/grafos_clustering.png')



#### Gráficos de silhouettes, con los mismos colores que los gráficos previos

sils = [silhouettes(dolph, particion) for particion in particiones]
sils_flat = [[s for cluster in sil_vals for s in cluster] for sil_vals in sils]
sils_prom = [sum(s)/len(s) for s in sils_flat]

fig, axes = plt.subplots(2, 4, figsize=(12, 6))
axes = np.ravel(axes)
for ax, metodo, sil_vals, colores, sprom in zip(axes, metodos, sils, colores_clusters, sils_prom):    
    graficar_silhouettes(sil_vals, colores=colores, ax=ax)
    titulo = metodo + r' ($\bar{s}$ = ' + '{:.2g})'.format(sprom)
    ax.set_title(titulo)
    # También guardamos el gráfico por separado
    fig_temp, ax_temp = plt.subplots()
    ax_temp.set_title(titulo)
    graficar_silhouettes(sil_vals, colores=colores, ax=ax_temp)
    fig_temp.tight_layout()
    fig_temp.savefig('Tp3/graficos/silhouette_corregido/Silhouettes {}.png'.format(metodo))
axes[7].axis('off')
fig.tight_layout()
fig.savefig('Tp3/graficos/silhouette_corregido/silhouettes_todos.png')