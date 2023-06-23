#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 30 03:16:02 2021

@author: david
"""


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
from scipy.spatial import Voronoi, voronoi_plot_2d, ConvexHull,KDTree
import math
from GeometriaTarea2 import FlightOperations


airports=pd.read_csv('airports_CO.dat',header=None, sep='\s\s+', engine='python',names=['latitude', 'longitude', 'altitude','Ciudad', 'departament', 'name'])
points=airports[['longitude','latitude']].to_numpy()
Borders=pd.read_csv('borders_CO.dat',header=None, sep='\s\s+', engine='python',names=['latitude', 'longitude'])



##Punto 1

newFlightOperations=FlightOperations(points,Borders,airports)


newFlightOperations.plot()

##Punto 2

newFlightOperations.Punto_2()


##Punto3

newFlightOperations.Punto_3()


##Punto 4

newFlightOperations.Punto_4()

#Punto 5

newFlightOperations.Punto_5()