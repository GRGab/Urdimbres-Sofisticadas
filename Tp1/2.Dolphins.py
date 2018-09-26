# -*- coding: utf-8 -*-
"""
Created on Fri Aug 31 18:39:12 2018

@author: Gabo
"""

import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from networkx.readwrite.gml import read_gml

from random import sample

from lectura import ldata
import sys
sys.path.append('./Tp1/')
from dolphins_funciones import (genero_a_color, particionar_por_genero,
                                crear_leyenda, contar_enlaces_internos,
                                contar_enlaces_entre_grupos,
                                contar_clases,p_value)
from modularidad import modularidad
from histograma import histograma
from graficar_multipartito import *
from collections import Counter


dolph = read_gml('Tp1/tc01_data/new_dolphins.gml')
genders = dict(ldata('Tp1/tc01_data/dolphinsGender.txt'))
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

contar_clases(dolph2, 'gender', ['f','m'])
# Hay 24 delfines hembra y 34 delfines macho.

#%%
n_simulaciones = int(1000)
enlaces_entre_grupos = np.zeros((n_simulaciones))
modularidades = np.zeros((n_simulaciones))
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
    modularidades[i] = modularidad(grafo_h0, 'gender')
    # Generar visualización para cada grafo (solo descomentar si
    # n_simulaciones es menor a 10!!!!))
#    fig, ax = plt.subplots()
#    colores = [genero_a_color(g) for g in nx.get_node_attributes(grafo_h0, "gender").values()]    
#    multi_spring_pos = position_multipartito_spring(grafo_h0, ['f', 'm'],
#                                                    'gender', dhorizontal=1.5)
#    nx.draw(grafo_h0, ax = ax, node_size = ns, node_color=colores,
#            pos=multi_spring_pos)
#    ax.set_title('Género aleatorizado, delfines sin información ignorados')
#    crear_leyenda(ax)
#%%
# Visualizar distribución de enlaces entre grupos bajo hipótesis nula
valor_real = contar_enlaces_entre_grupos(dolph2, 'gender')
fig, ax = histograma(enlaces_entre_grupos, bins=150, density=True,
                     titulo=r'Distribución de enlaces entre delfines de géneros distintos bajo $H_0$',
                     xlabel='# de enlaces')
ax.axvline(valor_real, color='deeppink',
           label='Valor real = {}'.format(valor_real))
ax.legend()
plt.show()
#%%
#Calculo del p-value
p_value(enlaces_entre_grupos)
#%%
# Visualizar distribución de modularidades
modularidad_real = modularidad(dolph2, 'gender')
fig, ax = histograma(modularidades, bins=15, density=True,
                     titulo=r'Distribución de modularidad bajo $H_0$',
                     xlabel='Modularidad')
ax.axvline(modularidad_real, color='deeppink',
           label='Valor real = {}'.format(valor_real))
     
#%%
#Intentos de Mati de hacer esto sin haber cursado estadistica. Riansenn.
#def p_value_1(datos, valor_real=52):     
#    a=[]
#    #Creo 2 listas una con los numeros de enlaces que aparecen y otra
#    # con cuantas veces se repiten cada numero enlace ("la altura de cada bin").
#    enlace = list(Counter(datos).keys())
#    cuantos = list(Counter(datos).values())
#    #Creo una lista comparador, cada componente sera la "altura del hist" si el
#    #valor medido se encuentra en mi histograma y sino un 1
#    comparador = []
#    for j in range (0,len(enlace)):
#        if enlace[j] == valor_real:
#            comparador.append(cuantos[j])
#        else:
#            comparador.append(1)
#    #Componente a componente, voy a comparar cada valor de "altura" con el comparador.
#    #Si es menor o igual al mismo, lo appendeo a una lista.
#    for i in range (0,len(enlace)):
#        if cuantos[i] <= comparador[i]:
#            a.append(cuantos[i])
#    #Finalemente, el P-value sera el largo de esta lista dividido el numero
#    #total de eventos
#    return len(a)/len(datos)
#
#def p_value_2(datos, bin1=0, bin2=1, nbins = 150, valor_real=52):
#    # use _ to assign the patches to a dummy variable since we don't need them
#    conteos, bins, _ = plt.hist(datos, nbins)    
#    # get the width of each bin
#    bin_width = bins[1] - bins[0]
#    cnorm = conteos / (np.sum(conteos) * bin_width)
#    conteos = cnorm    
#    # sum over number in each bin and mult by bin width, which can be factored out
#    integral = bin_width * sum(conteos[bin1:bin2])    
#    return integral


