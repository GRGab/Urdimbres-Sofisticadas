#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 13 02:04:34 2018

@author: gabo
"""

from lectura import ldata
import networkx as nx
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import random
import time
import traceback

import sys
sys.path.append('./Tp2/')
from agregar_esencialidad import agregar_esencialidad
from histograma import histograma
from funciones_de_desarme_ejc import (agregar_centralidad, 
                                        desarme_por_centralidad,
                                        desarme_por_centralidad_flow,
                                        desarme_por_centralidad_random)

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

nombres = {g_apms: "apms", g_lit: "lit", g_lit_reg: "lit_reg",
           g_y2h: "y2h"}
#%%

def log_traceback(ex, ex_traceback=None):
    if ex_traceback is None:
        ex_traceback = ex.__traceback__
    tb_lines = [ line.rstrip('\n') for line in
                 traceback.format_exception(ex.__class__, ex, ex_traceback)]
    for line in tb_lines:
        print(line)

#for g in [g_apms, g_lit, g_lit_reg, g_y2h]:
for g in [g_y2h]:
    for criterio in ['degree', 'eigen', 'sub', 'bet', 'flow', 'random']:
        if ((g is g_apms and criterio in ['degree', 'sub', 'flow']) or
            (g is g_lit and criterio is 'sub')):
            pass
        else:
            if criterio is 'flow':
                desarme = desarme_por_centralidad_flow
            elif criterio is 'random':
                desarme = desarme_por_centralidad_random
            else:
                desarme = desarme_por_centralidad
            
            print('Procesando g = ', nombres[g], ', crit = ', criterio)
            try:
                ti = time.time()
                cant_nodos_red, cant_nodos_cg = desarme(g, criterio=criterio)
                tf = time.time(); print(tf-ti, 'segundos')
                np.savez('curva_desarme_{}_{}.npz'.format(nombres[g], criterio),
                         cant_nodos_red=cant_nodos_red,
                         cant_nodos_cg=cant_nodos_cg)
                print('Guardado con éxito!')
            except Exception as ex:
                log_traceback(ex)
                print('La ejecución falló. Continuamos con el próximo caso.')
    print('Terminó el script.')