#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 28 18:34:25 2018

@author: matias
"""

from lectura import ldata
import networkx as nx
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd

from os.path import join as osjoin
import sys
sys.path.append('./Tp2/')
from agregar_esencialidad import agregar_esencialidad_dict
from histograma import histograma
from funciones_de_desarme_ejc import (agregar_centralidad, 
                                        desarme_por_centralidad,
                                        desarme_por_centralidad_flow,
                                        desarme_por_centralidad_random)
#%%
apms = ldata('Tp2/tc02Data/yeast_AP-MS.txt')

y2h = ldata('Tp2/tc02Data/yeast_Y2H.txt')

lit = ldata('Tp2/tc02Data/yeast_LIT.txt')

lit_r = ldata('Tp2/tc02Data/yeast_LIT_Reguly.txt')
lit_r = [fila[:2] for fila in lit_r[1:]]

g_apms = nx.Graph()
g_apms.add_edges_from(apms)

g_lit = nx.Graph()
g_lit.add_edges_from(lit)

g_lit_reg = nx.Graph()
g_lit_reg.add_edges_from(lit_r)

g_y2h = nx.Graph()
g_y2h.add_edges_from(y2h)

ess = ldata('Tp2/tc02Data/Essential_ORFs_paperHe.txt')
ess =  ess[2:-4]
ess = [fila[1] for fila in ess]
ess = np.unique(ess)

for g in [g_apms, g_lit, g_lit_reg, g_y2h]:
    agregar_esencialidad_dict(g, ess)
#%%
#    Ejercico c
#Punto I

nombres = {g_apms: "apms", g_lit: "lit", g_lit_reg: "lit_reg",
           g_y2h: "y2h"}
curvas = {}
for g in [g_apms, g_lit, g_lit_reg, g_y2h]:
    for criterio in ['degree', 'eigen', 'sub', 'bet', 'flow', 'random']:
        filename = 'curva_desarme_{}_{}.npz'.format(nombres[g], criterio)
        data = np.load(osjoin('Tp2', 'tc02Data', filename))
        curvas['{}_{}'.format(nombres[g], criterio)] = [data['cant_nodos_red'],
                                                        data['cant_nodos_cg']]
        
# Además de las curvas, queremos el puntito que nos marca a dónde va a parar
# el tamaño de la cg cuando eliminamos de un saque todas las esenciales

def desarme_esenciales_saco_ess(g,ess):
    '''
    Toma una red y una lista de nodos esenciales. Elimina los nodos esenciales
    y calcula el tamaño relativo de la componente gigante en la nueva red.
    '''
   #Parametros para la normalizacion
    cg_original = max(nx.connected_component_subgraphs(g), key=len).order()
    #Copio la red para hacerla bolsa
    G = g.copy()
    G.remove_nodes_from(ess)
    #calculo la componente gigante y su tamaño
    cg = max(nx.connected_component_subgraphs(G), key=len)
    maxcomp = cg.order()
    a = maxcomp/cg_original #Lo que da sacando esenciales
    return a
#%% Graficar para cada red por separado
        
fontsize = 18
ticksize = 16
        
for g in [g_apms, g_lit, g_lit_reg, g_y2h]:
    with plt.style.context(('seaborn')):
        fig, ax = plt.subplots(figsize=(10,8))
        
    # Dibujamos las curvas
    for criterio in ['degree', 'eigen', 'sub', 'bet', 'flow', 'random']:
        xs = curvas['{}_{}'.format(nombres[g], criterio)][0]
        ys = curvas['{}_{}'.format(nombres[g], criterio)][1]
        if criterio is 'random': # Promedio las 100 historias
            xs = np.average(xs, axis=0)
            ys = np.average(ys, axis=0)
            ax.plot(xs, ys, '-', label=criterio)
        else:
            ax.plot(xs, ys, '.-', label=criterio)
            
    ### Dibujamos puntito de remoción de esenciales
    # Cuento nodos esenciales en esta red particular
    num_esenciales = sum(dict(g.nodes(data='esencialidad')).values())
    # Calculo fracción
    fracnodos_esenciales = num_esenciales / g.order() # Esto habría que cambiarlo si rehiciéramos las curvas con las cg únicamente
    # Calculo tamaño de la cg al eliminarlos
    cg_sinesenciales = desarme_esenciales_saco_ess(g, ess)
    # Dibujo
    ax.plot(fracnodos_esenciales, cg_sinesenciales, '*', label='Esenciales',
            ms=10)
    
    # Emprolijamos
    ax.tick_params(labelsize=ticksize)
    ax.set_xlabel('Fracción de nodos removidos',
              fontsize=fontsize)
    ax.set_ylabel('Tamaño relativo de la componente gigante',
              fontsize=fontsize)
    ax.legend(fontsize=fontsize)
    ax.set_title(nombres[g], fontsize=fontsize)
    ax.set_xlim([-0.01, 0.45])
    fig.tight_layout()
    fig.savefig('Tp2/Ej c/curvas_desarme_{}.png'.format(nombres[g]))