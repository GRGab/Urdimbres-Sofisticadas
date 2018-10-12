from __future__ import division
import numpy as np
import random
import copy
import time

def random_edges(edges):
    new_edges = copy.deepcopy(edges)
    sequence = [edges[i][1] for i in range(len(edges))]
    np.random.shuffle(sequence)
    for i in range(len(edges)):
        new_edges[i][1] = sequence[i]
    return new_edges

edges = [[1, 2], [1, 3], [2, 4], [2, 3], [3, 4]]

d_edges = random_edges(edges)


t = time.time()
d_edges = random_edges(y2h)
elapsed = time.time() - t
print(elapsed)

def count_interactions(G, ess_nodes):
    sub_G = G.subgraph(ess_nodes)
    return sub_G.size()

t = time.time()
datos_redes = [apms, lit, y2h, lit_r]
redes = [g_apms, g_lit, g_y2h, g_lit_reg]
m_original = []
m_mean = []
for h in range(len(redes)):
    m = []
    red_size = redes[h].size()
    for i in range(10000):
        d_edges = random_edges(datos_redes[h])
        g = nx.Graph()
        g.add_edges_from(d_edges)
        m.append(count_interactions(g, ess))

    m_original.append(count_interactions(redes[h], ess))
    m_mean.append(np.mean(m))
elapsed = time.time() - t
print(elapsed)    
    
m_original = count_interactions(g_apms, ess)
fig, ax = histograma(m, bins=100, density=True, errbars=False, 
                     titulo=r'Distribucion de enlaces escenciales para recableos de la original',
                     xlabel='Numero de enlaces escenciales')
ax.axvline(m_original, color='deeppink',
           label='Valor real = {}'.format(m_original))
ax.legend()
plt.show()

#plt.savefig('histograma_he')
m = np.loadtxt('/home/tomas/Desktop/Redes complejas/Urdimbres-Sofisticadas/datos_para_histograma_He_apms')

m_mean = np.mean(m)

alpha_apms = (m_original - m_mean)/g_apms.size()

m = m[51]
G = g_apms.copy()
g, values = agregar_esencialidad(g_apms, ess)
len_nodes_ess_original =  np.sum(values)
edges_new_ess = random.sample(G.edges(), int(m_original-m))
nodos_ess_temp = [[edges_new_ess[i][1] for i in range(len(edges_new_ess))], [edges_new_ess[i][0] for i in range(len(edges_new_ess))]]
nodos_ess = [item for sublist in nodos_ess_temp for item in sublist]
nodos = np.unique(nodos_ess)

new_graph_ess = G.subgraph(edges_new_ess)
nodos_ess_temp = [[edges_new_ess[i][1] for i in range(len(edges_new_ess))], [edges_new_ess[i][0] for i in range(len(edges_new_ess))]]
nodos_ess = [item for sublist in nodos_ess_temp for item in sublist]
G.remove_nodes_from(nodos_ess)
new_edges_ess_random = G.subgraph(random.sample(G.edges(), int(m)))

