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
d_edges = random_edges(apms)
elapsed = time.time() - t
print(elapsed)

def count_interactions(G, ess_nodes):
    sub_G = G.subgraph(ess_nodes)
    return sub_G.size()

t = time.time()
m = []
for i in range(10000):
    d_edges = random_edges(apms)
    g = nx.Graph()
    g.add_edges_from(d_edges)
    m.append(count_interactions(g, ess))
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

plt.savefig('histograma_he')
np.savetxt('datos_para_hisograma_He_apms', m)