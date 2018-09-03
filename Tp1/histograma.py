# -*- coding: utf-8 -*-
"""
Created on Sun Sep  2 21:38:04 2018

@author: Gabo
"""

import numpy as np
import matplotlib.pyplot as plt

#Estilos disponibles para pyplot:
#['bmh', 'classic', 'dark_background', 'fast', 'fivethirtyeight', 'ggplot',
# 'grayscale', 'seaborn-bright', 'seaborn-colorblind', 'seaborn-dark-palette',
# 'seaborn-dark', 'seaborn-darkgrid', 'seaborn-deep', 'seaborn-muted',
# 'seaborn-notebook', 'seaborn-paper', 'seaborn-pastel', 'seaborn-poster',
# 'seaborn-talk', 'seaborn-ticks', 'seaborn-white', 'seaborn-whitegrid',
# 'seaborn', 'Solarize_Light2', '_classic_test']

def histograma(valores_a_binear, bins=None, titulo_hist=None, magnitud_x=None,
               density=False):
    with plt.style.context(('seaborn')):
            fig, ax = plt.subplots()
            
    conteos, bordes_bines, _ = plt.hist(valores_a_binear, bins=bins,
                                        density=density)

    # Graficar errores
    error = np.sqrt(conteos)
    # Opcional: normalizar el error si se normalizan los conteos:
    if density==True:
        error = np.sqrt(conteos)/np.sqrt(len(valores_a_binear))
    # OJO: esta fórmula es incorrecta si el ancho de los bines no es constante
    centros_bines = 0.5 * (bordes_bines[:-1] + bordes_bines[1:])
    plt.errorbar(centros_bines, conteos,
                 yerr=error,
                 fmt='none',
                 capsize=0,
                 ecolor='k')
    
    if titulo_hist != None:
        ax.set_title(titulo_hist, fontsize=16)
    if magnitud_x != None:
        ax.set_xlabel(magnitud_x, fontsize=14)
    ylabel = '# de eventos' if density==False else '# de eventos normalizado'
    ax.set_ylabel(ylabel, fontsize=14)
    num_bines = len(bordes_bines) - 1
    anotacion = ('$N = $' + str(len(valores_a_binear))+ '\n' +
                 r'$N_{bines}$ = ' + str(num_bines))
    ax.annotate(anotacion,
                (.8, .8), xycoords='axes fraction',
                backgroundcolor='w', fontsize=14)
    fig.tight_layout()
    
    # Guardar imagen
    #plt.savefig(os.path.join(graficos_path, titulo_hist) + '.png')
    plt.show()
    return fig, ax