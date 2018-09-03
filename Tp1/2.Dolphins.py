# -*- coding: utf-8 -*-
"""
Created on Fri Aug 31 18:39:12 2018

@author: Gabo
"""

import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from networkx.readwrite.gml import read_gml
from graficar_multipartito import position_multipartito, position_multipartito_spring
from lectura import ldata
from random import sample
from dolphins_funciones import (genero_a_color, particionar_por_genero,
                                crear_leyenda, contar_enlaces_internos,
                                contar_enlaces_entre_grupos)
from histograma import histograma



dolph = read_gml('tc01_data/new_dolphins.gml')
genders = dict(ldata('tc01_data/dolphinsGender.txt'))
#%%

# Agrego los sexos a los dicts de cada delfín
for nodo, dict_nodo in dict(dolph.nodes).items():
    dict_nodo['gender'] = genders[nodo] # agrego el sexo del delfín a su dict
#    print('Key = {}, Value = {}'.format(nodo, dict_nodo)) # para chequear que anda


particiones = particionar_por_genero(dolph)
colores = [genero_a_color(g) for g in nx.get_node_attributes(dolph, "gender").values()]    

# Hay 34 delfines macho, 24 delfines hembra y 4 delfines sin información

#%%
# Graficamos de 9 maneras diferentes

fig, axes = plt.subplots(3,3)
axes = axes.flatten()
ns = 25 # Node size
# Posicionamiento multipartito (es engañoso: no se ven los links homofílicos)
multipartite_pos = position_multipartito(dolph, ['f', 'NA', 'm'], 'gender')
nx.draw(dolph, ax = axes[0], node_size = ns, node_color=colores,
        pos=multipartite_pos)
# Posicionamiento en un círculo
nx.draw_circular(dolph, ax = axes[1], node_size = ns, node_color=colores)
# Posicionamiento en círculos concéntricos
shell_pos = nx.drawing.layout.shell_layout(dolph, particiones)
nx.draw(dolph, ax = axes[2], node_size = ns, node_color=colores,
        pos=shell_pos)
# Posicionamiento al azar
nx.draw_random(dolph, ax = axes[3], node_size = ns, node_color=colores)
# Posicionamiento espectral
nx.draw_spectral(dolph, ax = axes[4],
                 node_size = ns, node_color=colores)
# Posicionamiento por resortes
nx.draw_spring(dolph, ax = axes[5],
                 node_size = ns, with_labels=False, node_color=colores)
# Posicionamiento multipartito al azar. Posiciono al azar y
# luego desplazo lateralmente según género
multi_random_pos = position_multipartito_random(dolph,
                                                ['f', 'm', 'NA'], 'gender')
nx.draw(dolph, ax = axes[6], node_size = ns, node_color=colores,
        pos=multi_random_pos)
# Posicionamiento multipartito espectral. Posiciono por espectro y
# luego desplazo lateralmente según género
multi_spectral_pos = position_multipartito_spectral(dolph, ['f', 'm', 'NA'],
                                                    'gender', dhorizontal=0.5)
nx.draw(dolph, ax = axes[7], node_size = ns, node_color=colores,
        pos=multi_spectral_pos)
# Posicionamiento multipartito por resortes. Posiciono por resortes y
# luego desplazo lateralmente según género
multi_spring_pos = position_multipartito_spring(dolph, ['f', 'm', 'NA'],
                                                'gender', dhorizontal=2)
nx.draw(dolph, ax = axes[8], node_size = ns, node_color=colores,
        pos=multi_spring_pos)


#%%
# 3 Gráficos últimos gráficos (onda bipartitos) de la figura anterior:
ns = 25

############ Random con desplazamiento lateral por género #
fig, ax = plt.subplots()
multi_random_pos = position_multipartito_random(dolph,
                                                ['f', 'm', 'NA'], 'gender')
nx.draw(dolph, ax = ax, node_size = ns, node_color=colores,
        pos=multi_random_pos)
ax.set_title('Layout aleatorio con desplazamiento lateral por género')
crear_leyenda(ax)
###########################################################

########## Espectral con desplazamiento lateral por género #
fig, ax = plt.subplots()
multi_spectral_pos = position_multipartito_spectral(dolph, ['f', 'm', 'NA'],
                                                    'gender', dhorizontal=0.5)
nx.draw(dolph, ax = ax, node_size = ns, node_color=colores,
        pos=multi_spectral_pos)
ax.set_title('Layout por espectro con desplazamiento lateral por género')
crear_leyenda(ax)
###########################################################

########## Resortes con desplazamiento lateral por género #
fig, ax = plt.subplots()
multi_spring_pos = position_multipartito_spring(dolph, ['f', 'm', 'NA'],
                                                'gender', dhorizontal=1.5)
nx.draw(dolph, ax = ax, node_size = ns, node_color=colores,
        pos=multi_spring_pos)
ax.set_title('Layout por resortes con desplazamiento lateral por género')
crear_leyenda(ax)
###########################################################

#%%

# PUNTO B: Cuantificar homofilia

# B) I) generamos muchas redes a partir de nuestra red real redistribuyendo
# los géneros de manera aleatoria. Respetamos las cantidades de delfines de
# cada género. Dejamos de lado a los delfines cuyo género no conocemos.

delfines_con_info = [d for d in dolph.nodes() if d not in particiones[1]]
dolph2 = nx.subgraph(dolph, delfines_con_info).copy()
# De ahora en adelante trabajamos con dolph2 en vez de dolph
# Hay 24 delfines hembra y 34 delfines macho.


n_simulaciones = int(1e4)
enlaces_entre_grupos = np.zeros((n_simulaciones))
grafo_h0 = dolph2.copy()
# Vamos a ir modificando este grafo "in place" (no lo clonamos n veces)

for i in range(n_simulaciones):
    # Mezclamos la lista de nombres de delfines.
    # A los primeros 24 delfines les reasignamos género hembra
    # El resto van a ser macho.
    hembras = sample(list(grafo_h0.nodes()), 24)
    for nombre in grafo_h0.nodes():
        if nombre in hembras:
            grafo_h0.nodes()[nombre]['gender'] = 'f'
        else:
            grafo_h0.nodes()[nombre]['gender'] = 'm'
    enlaces_entre_grupos[i] = contar_enlaces_entre_grupos(grafo_h0, 'gender')
    
    # Generar visualización para cada grafo (solo descomentar si
    # n_simulaciones es menor a 10!!!!))
    # fig, ax = plt.subplots()
    # colores = [genero_a_color(g) for g in nx.get_node_attributes(grafo_h0, "gender").values()]    
    # multi_spring_pos = position_multipartito_spring(grafo_h0, ['f', 'm'],
    #                                                 'gender', dhorizontal=1.5)
    # nx.draw(grafo_h0, ax = ax, node_size = ns, node_color=colores,
    #         pos=multi_spring_pos)
    # ax.set_title('Género aleatorizado, delfines sin información ignorados')
    # crear_leyenda(ax)
#%%
# Visualizar distribución de enlaces entre grupos bajo hipótesis nula
fig, ax = histograma(enlaces_entre_grupos, bins=15, density=True,
                     titulo_hist=r'Distribución de enlaces entre delfines de géneros distintos bajo $H_0$',
                     magnitud_x='# de enlaces')
valor_real = contar_enlaces_entre_grupos(dolph2, 'gender')
ax.axvline(valor_real, color='deeppink',
           label='Valor real = {}'.format(valor_real))
ax.legend()