import numpy as np
import networkx as nx
from networkx.readwrite.gml import read_gml

import sys
sys.path.append('./Tp3/')
from funciones_tp3 import calcular_particion, comunidad_a_color
import matplotlib.pyplot as plt

#%%

dolph = read_gml('Tp3/dolphins.gml')    
lista = ["infomap","label_prop", "fastgreedy", "eigenvector", "louvain"
        , "edge_betweenness", "walktrap"]


npzfile = np.load('Tp3/tc03Data/Ej_b_particiones_tomi.npz')
rewire = npzfile['salida']
original = npzfile['salida_grafo_original']
## EJ 1 A ... Visualizacion de las distintas particiones.

colors = []
for metodo in lista:
    nodes = calcular_particion(dolph, method = metodo)
    colors.append(comunidad_a_color(dolph, nodes))
#[ax1, ax2, ax3, ax4], [ax5, ax6, ax7, ax8]

ns = 35
#nx.draw(dolph, node_color = colors[0])
pos = nx.spring_layout(dolph)

fig, ([ax1, ax2, ax3, ax4], [ax5, ax6, ax7, ax8]) = plt.subplots(2, 4)

plt.sca(ax1)
ax1.set_title('Infomap')
nx.draw(dolph, node_size = ns, node_color=colors[0], pos = pos)

plt.sca(ax2)
ax2.set_title('Label prop')
nx.draw(dolph, node_size = ns, node_color=colors[1], pos = pos)

plt.sca(ax3)
ax3.set_title('Fastgreedy')
nx.draw(dolph, node_size = ns, node_color=colors[2], pos = pos)

plt.sca(ax4)
ax4.set_title('Eigenvectors')
nx.draw(dolph, node_size = ns, node_color=colors[3], pos = pos)

plt.sca(ax5)
ax5.set_title('Louuvain')
nx.draw(dolph, node_size = ns, node_color=colors[4], pos = pos)

plt.sca(ax6)
ax6.set_title('Edge betweenness')
nx.draw(dolph, node_size = ns, node_color=colors[5], pos = pos)

plt.sca(ax7)
ax7.set_title('Walk trap')
nx.draw(dolph, node_size = ns, node_color=colors[6], pos = pos)

empty = nx.Graph()
plt.sca(ax8)
nx.draw(empty)