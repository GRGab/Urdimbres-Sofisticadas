# -*- coding: utf-8 -*-
"""
Created on Fri Aug 31 18:39:12 2018

@author: Gabo
"""

import networkx as nx
import matplotlib.pyplot as plt
from networkx.readwrite.gml import read_gml
from lectura import ldata

dolph = read_gml('tc01_data/new_dolphins.gml')
genders = dict(ldata('tc01_data/dolphinsGender.txt'))
#%%

# Agrego los sexos a los dicts de cada delfín
for key, value in dict(dolph.nodes).items():
    dict_nodo = dolph.nodes[key] # agarro el dict de un delfín
    dict_nodo['gender'] = genders[key] # agrego el sexo del delfín a su dict
#    print('Key = {}, Value = {}'.format(key, value)) # para chequear que anda

def genero_a_color(gender):
    if gender=='m':
        return 'green'
    elif gender=='f':
        return 'dodgerblue'
    else:
        return 'black'
    
def particionar_por_genero(G, orden={'f':2, 'm':0, 'NA':1}):
    particiones = [[], [], []]
    for key in dict(G.nodes).keys():
        gender = G.nodes[key]['gender']
        if gender=='f':
            particiones[orden[gender]].append(key)
        elif gender=='m':
            particiones[orden[gender]].append(key)
        else:
            particiones[orden[gender]].append(key)
    return particiones

colores = [genero_a_color(g) for g in nx.get_node_attributes(dolph, "gender").values()]    

fig, axes = plt.subplots(2,2)
axes = axes.flatten()
nx.draw_circular(dolph, ax = axes[0],
                 node_size = 50, with_labels=False, node_color=colores)
shell_pos = nx.drawing.layout.shell_layout(dolph,
                                           particionar_por_genero(dolph))
nx.draw(dolph, ax = axes[1], node_size = 50, node_color=colores,
        pos=shell_pos)
nx.draw_spectral(dolph, ax = axes[2],
                 node_size = 50, node_color=colores)
nx.draw_spring(dolph, ax = axes[3],
                 node_size = 50, with_labels=False, node_color=colores)

