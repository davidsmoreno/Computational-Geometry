#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  7 09:38:50 2021

@author: david
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay
import matplotlib.tri as mtri
from scipy.spatial import Delaunay

from scipy.interpolate import LinearNDInterpolator


import pandas as pd
from GeometriaTarea1 import TIN

df=pd.read_table("pts1000c.dat",sep=' ')
df=df.to_numpy()

puntos=df[:,[0,1]]
z = df[:,2]
num_puntos=len(puntos)

#Creamos la clase

newTIN=TIN(num_puntos,puntos,z)

#Ploteamos la triangulación en 2d con la altua pintada de color rojo

newTIN.graph_punto1_2d()

newTIN.graph_punto1_3d()

##Punto 2 Calcula la altura dada la interpolacion del triangulo en el que está inscrito
punto=[1,1]

newTIN.Interpolacion_punto2(punto)

##Grafica punto
newTIN.gratafica_Interpolacion_punto2(punto)


##Punto 3

Vertices_Cuadrilatero=newTIN.Largest_area(punto)
print(Vertices_Cuadrilatero)

#Grafica
newTIN.grafica_Largest_area(Vertices_Cuadrilatero)



#Punto 4 
L=newTIN.max_min_angulos()

L2=newTIN.angulos_tringulacion()


#Punto 5
S_T=newTIN.SpanningTree()
print(S_T,"Spanning Tree")




