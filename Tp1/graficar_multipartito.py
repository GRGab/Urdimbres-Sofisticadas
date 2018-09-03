# -*- coding: utf-8 -*-
"""
Created on Sun Sep  2 19:05:25 2018

@author: Gabo
"""

import networkx as nx
import numpy as np


def position_multipartito(g, clases, nombre_clase):
  """g es un objeto Graph de networkX, cuyos nodos deben tener un atributo
  llamado nombre_clase. clases es una lista con los valores posibles que
  puede adoptar dicho atributo (que corresponderán cada uno a una columna).
  Devuelve un diccionario cuyos keys son los nombres de los nodos y cuyos
  valores son las posiciones para graficar.
  Ejemplo con datos de los delfines:
      clases == ['f', 'm', 'NA']
      nombre_clase == 'gender'
      
  Modificado a partir de código de James A. Foster en
  https://stackoverflow.com/questions/35472402/how-do-display-bipartite-graphs-with-python-networkx-package
  """
  xPos = {}
  yPos = {}
  for indice, clase in enumerate(clases):
     xPos[clase] = indice
     yPos[clase] = 0

  pos = {}
  for node, attrDict in g.nodes(data=True): # Piola sintaxis lo de data=True
     clase_nodo = attrDict[nombre_clase]
     # print('Nodo: {}\t{}: {}'.format(node, nombre_clase, clase_nodo))     
     # print('\t(x,y): ({},{})'.format(xPos[clase_nodo], yPos[clase_nodo]))
     pos[node] = (xPos[clase_nodo], yPos[clase_nodo])
     yPos[clase_nodo] += 2

  return pos

def position_multipartito_random(g, clases, nombre_clase, dhorizontal=2):
  """g es un objeto Graph de networkX, cuyos nodos deben tener un atributo
  llamado nombre_clase. clases es una lista con los valores posibles que
  puede adoptar dicho atributo (que corresponderán cada uno a una columna).
  dhorizontal es el espaciado entre clases diferentes de nodos.
  Devuelve un diccionario cuyos keys son los nombres de los nodos y cuyos
  valores son las posiciones para graficar.
  Ejemplo con datos de los delfines:
      clases == ['f', 'm', 'NA']
      nombre_clase == 'gender'
  """
  
  pos = nx.random_layout(g)
  for node, attrDict in g.nodes(data=True):
      clase_nodo = attrDict[nombre_clase]
      for i, clase in enumerate(clases):
          if clase_nodo==clase:
              adicion = dhorizontal*i
      posx_old = pos[node][0]
      posy_old = pos[node][1]
      pos[node] = np.array([posx_old + adicion, posy_old])
  return pos

def position_multipartito_spectral(g, clases, nombre_clase, dhorizontal=2):
  """g es un objeto Graph de networkX, cuyos nodos deben tener un atributo
  llamado nombre_clase. clases es una lista con los valores posibles que
  puede adoptar dicho atributo (que corresponderán cada uno a una columna).
  dhorizontal es el espaciado entre clases diferentes de nodos.
  Devuelve un diccionario cuyos keys son los nombres de los nodos y cuyos
  valores son las posiciones para graficar.
  Ejemplo con datos de los delfines:
      clases == ['f', 'm', 'NA']
      nombre_clase == 'gender'
  """
  
  pos = nx.spectral_layout(g)
  for node, attrDict in g.nodes(data=True):
      clase_nodo = attrDict[nombre_clase]
      for i, clase in enumerate(clases):
          if clase_nodo==clase:
              adicion = dhorizontal*i
      posx_old = pos[node][0]
      posy_old = pos[node][1]
      pos[node] = np.array([posx_old + adicion, posy_old])
  return pos

def position_multipartito_spring(g, clases, nombre_clase, dhorizontal=1):
  """g es un objeto Graph de networkX, cuyos nodos deben tener un atributo
  llamado nombre_clase. clases es una lista con los valores posibles que
  puede adoptar dicho atributo (que corresponderán cada uno a una columna).
  dhorizontal es el espaciado entre clases diferentes de nodos.
  Devuelve un diccionario cuyos keys son los nombres de los nodos y cuyos
  valores son las posiciones para graficar.
  Ejemplo con datos de los delfines:
      clases == ['f', 'm', 'NA']
      nombre_clase == 'gender'
  """
  
  pos = nx.spring_layout(g)
  for node, attrDict in g.nodes(data=True):
      clase_nodo = attrDict[nombre_clase]
      for i, clase in enumerate(clases):
          if clase_nodo==clase:
              adicion = dhorizontal*i
      posx_old = pos[node][0]
      posy_old = pos[node][1]
      pos[node] = np.array([posx_old + adicion, posy_old])
  return pos

def position_multipartito_kk(g, clases, nombre_clase, dhorizontal=1):
  """g es un objeto Graph de networkX, cuyos nodos deben tener un atributo
  llamado nombre_clase. clases es una lista con los valores posibles que
  puede adoptar dicho atributo (que corresponderán cada uno a una columna).
  dhorizontal es el espaciado entre clases diferentes de nodos.
  Devuelve un diccionario cuyos keys son los nombres de los nodos y cuyos
  valores son las posiciones para graficar.
  Ejemplo con datos de los delfines:
      clases == ['f', 'm', 'NA']
      nombre_clase == 'gender'
  """
  
  pos = nx.kamada_kawai_layout(g)
  for node, attrDict in g.nodes(data=True):
      clase_nodo = attrDict[nombre_clase]
      for i, clase in enumerate(clases):
          if clase_nodo==clase:
              adicion = dhorizontal*i
      posx_old = pos[node][0]
      posy_old = pos[node][1]
      pos[node] = np.array([posx_old + adicion, posy_old])
      
  return pos