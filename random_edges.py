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
d_edges = random_edges(apms)
g = nx.Graph()
g.add_edges_from(d_edges)
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
    red_size = g_apms.size()
    for i in range(1000):
        d_edges = random_edges(apms)
        g = nx.Graph()
        g.add_edges_from(d_edges)
        m.append(count_interactions(g, ess))

    m_original.append(count_interactions(g_lit_reg, ess))
    m_mean.append(np.mean(m))
elapsed = time.time() - t
print(elapsed)    
    
m_mean = np.mean(m)
m_original = count_interactions(g_apms, ess)
fig, ax = histograma(m, bins=105, density=True, errbars=True, 
                     titulo=r'Distribucion de enlaces escenciales para recableos de la original',
                     xlabel='Numero de enlaces escenciales')
ax.axvline(m_original, color='deeppink',
           label='Valor real = {}'.format(m_original))
ax.legend()
plt.show()

#plt.savefig('histograma_he')
m = np.loadtxt('/home/tomas/Desktop/Redes complejas/Urdimbres-Sofisticadas/Tp2/datos_para_histograma_He_apms')


alpha_apms = (m_original - m_mean)/g_apms.size()



G = g_lit.copy()
g, values = agregar_esencialidad(g_apms, ess)
len_nodes_ess_original = np.sum(values)
num_nodes = len(G.nodes())
beta = []
for i in range(len(m)):
    edges_new_ess = random.sample(G.edges(), int(m_original-m[i]))
    nodos_ess_temp = [[edges_new_ess[i][1] for i in range(len(edges_new_ess))], [edges_new_ess[i][0] for i in range(len(edges_new_ess))]]
    nodos_ess = [item for sublist in nodos_ess_temp for item in sublist]
    nodos = np.unique(nodos_ess)
#    nodos_relleno = copy.deepcopy(nodos)
#    nodos_relleno = list(nodos_relleno)
    nodos_relleno = []
    i = len(nodos)
    while i <= int(len_nodes_ess_original):
        k = random.choice(list(g_apms.nodes()))
        if  k not in nodos:
            i = i + 1
            nodos_relleno.append(k)
        else: 
            i = i
            nodos_relleno.append(k)
    beta.append(len(nodos_relleno)/float(num_nodes))

plt.hist(beta, bins=50)

beta_mean = np.mean(beta) #Comparar con beta_lit_recta! Dan re parecidos


out.beta[0]* np.exp(out.beta[0]) * out.sd_beta[0]
out.beta[1]* np.exp(out.beta[1]) * out.sd_beta[1]






