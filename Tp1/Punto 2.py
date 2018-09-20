#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 31 13:51:21 2018

@author: matias
"""

import networkx as nx
import numpy as np
from matplotlib import pyplot as plt
from lectura import ldata

dolphin = nx.read_gml('tc01_data/new_dolphins.gml')
genders = ldata('tc01_data/dolphinsGender.txt')

color = []

for i in range (0,len(genders)):
    dolphin.add_node(genders[i][0], gender = genders[i][1])
    if genders[i][1] == 'm':
        color.append('blue')
    elif genders[i][1]== 'f':
        color.append('red')
    else:
        color.append('green')
#
#fig, axes= plt.subplots(2,2)
#axes = axes.flatten()
#nx.draw_shell(dolphin, with_labels=True, node_color=color, node_size = 10
#              , ax=axes[0])

N = []
for i in range(len(genders)):
    N.append(genders[i][1])

shells1 = [genders[i][0] for i in range(len(N)) if genders[i][1] == 'NA']
shells2 = [genders[i][0] for i in range(len(N)) if genders[i][1] == 'm']
shells3 = [genders[i][0] for i in range(len(N)) if genders[i][1] == 'f']

shells= [shells1, shells2, shells3]
pos = nx.drawing.layout.shell_layout(dolphin, shells)

nx.draw(dolphin, with_labels=True, node_color=color, node_size = 20, pos=pos)
