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

def histograma(valores_a_binear, bins='auto', titulo=None, magnitud_x=None,
               density=False, logbins=False, logx=False, logy=False, ax=None):
    """Si logbins=True, el parámetro bins debe ser una tupla que contenga
    los dos extremos del intervalo que se desea binear, y el número de bines
    logarítmicos que se desea usar, en ese orden."""
    
    if ax is None:
        with plt.style.context(('seaborn')):
            fig, ax = plt.subplots()
    else:
        fig = ax.get_figure()
    
    if logbins:
        assert isinstance(bins, tuple), 'bins debe ser una tupla de 3 elementos'
        logstart = np.log10(bins[0])
        logstop = np.log10(bins[1])
        nbins = bins[2]
        bins = np.logspace(logstart, logstop, num=nbins)
        bins = np.concatenate(([0], bins))
        
    conteos, bordes_bines = np.histogram(valores_a_binear, bins=bins)
    w = np.diff(bordes_bines) # Anchos de los bines
    
    # Normalizar, si es necesario, y asignar errores
    if density == True:
        cnorm = conteos / (np.sum(conteos) * w)
        enorm = np.sqrt(conteos) / (np.sum(conteos) * w)
        conteos, errores = cnorm, enorm
    else:
        errores = np.sqrt(conteos)
    
    # Graficar
    ax.bar(bordes_bines[:-1], conteos, width=w, yerr=errores, align='edge',
           color='dodgerblue', capsize=0)
    
#    centros_bines = 0.5 * (bordes_bines[:-1] + bordes_bines[1:])
#    plt.errorbar(centros_bines, conteos,
#                 yerr=errores,
#                 fmt='none',
#                 capsize=0,
#                 ecolor='k')
    
    if titulo != None:
        ax.set_title(titulo, fontsize=16)
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
    
    plt.show()
    return fig, ax

def binplot(valores, imin=None, imax=None, titulo=None,
            errorbars=True, ax=None):
    """Gráfico de barras. Recibe un iterable con
    los valores a graficar, de forma tal que la altura del bin i es igual a
    valores[i]. imin e imax son los índices entre los cuales graficar. Si no se
    especifican, se grafica la región en la cual las alturas de los bines
    son disintas de cero.
        Útil para realizar histogramas de cantidades discretas.
    
    PENDIENTES:
        - Agregar opción de normalización
        - Agregar opción de graficar los errores de manera distinta (simétrica)
        en el caso de que log == True
        - Agregar opción `multiplo` que permita que en el eje x aparezcan
        múltiplos de un cierto valor (indicando las unidades en el rótulo)."""
    if imin is None:
        imin = np.where(valores != 0)[0][0]
        # Si puedo, corro imin 1 bin a la izquierda para que quede más lindo:
        if imin != 0:
            imin -= 1
    if imax is None:
        imax = np.where(valores != 0)[0][-1] + 1
        # El +1 es porque valores[imax] no será graficado.
        # Si puedo, corro imax 1 bin a la derecha para que quede más lindo:
        if imax != len(valores):
            imax += 1
    if ax is None:
        with plt.style.context(('seaborn')):
            fig, ax = plt.subplots()
    else:
        fig = ax.get_figure()
    xs, hs = range(imin, imax), valores[imin:imax]
    errs = np.sqrt(hs) if errorbars == True else None
    ax.bar(xs, hs, width=1, yerr=errs, align='edge',
           color='dodgerblue', capsize=0)
    if log == True:
        ax.set_yscale('log')
    ax.set_xticks(np.arange(imin, imax))
    if titulo is not None:
        ax.set_title(titulo, fontsize=16)
    fig.tight_layout()
    
def hist_discreto(xs, imin=None, imax=None, titulo=None,
                  log=False, errorbars=True, ax=None):
    """Realiza un histograma de una variable discreta, a partir
    de los samples contenidos en `xs`.
    
    Es un envoltorio para la función binplot.
    """
    hs = np.bincount(xs)
    binplot(hs, imin=imin, imax=imax, titulo=titulo, log=log,
            errorbars=errorbars, ax=ax)

    
if __name__ == '__main__':
    from scipy.stats import expon
    xs = expon(0.5).rvs(1000)
#    fig, ax = histograma(xs, titulo='El título', magnitud_x='El eje x')
#    histograma(xs, density=True)
    histograma(xs, logbins=True, bins=(1,10,20))
