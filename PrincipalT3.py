#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 25 23:04:30 2021

@author: david
"""

from GeometriaTarea3 import FlightOperations
import numpy as np
import random 
from scipy import spatial
from scipy.spatial import KDTree
import matplotlib.pyplot as plt 
from matplotlib.patches import Rectangle




##Semilla
random.seed(10)


Data_y = np.random.uniform(-1, 1, 9)
Data_x=np.random.uniform(-1, 1, 9)

##Creamos la nueva clase con las coordenadas (x,y) de los barcos
New_FlightOperations=FlightOperations(Data_x,Data_y)


##Punto 1

point=(0.2,0.5) ##Vessel a elegir en el intervalo x e (-1,1), yx (-1,1)

New_FlightOperations.Punto_1(point)

##Punto 2


radio=1
New_FlightOperations.Punto_2(radio)


##Punto 3
New_FlightOperations.punto_3()


##Punto4
New_FlightOperations.Plot_Kdtree()