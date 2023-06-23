#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 01:08:39 2021

@author: david
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt 
from scipy.spatial import Voronoi, voronoi_plot_2d, ConvexHull,KDTree
import math



class FlightOperations:
    
        
    def __init__(self,points,borders,airports):
        
        self.points=points
        self.borders=borders
        self.n_points=len(points)
        self.vor=Voronoi(self.points)
        self.airports=airports
        


    def plot(self):
        Borders=self.borders.to_numpy()
        #Grafica de Voronoi con sus Fronteras
        fig=voronoi_plot_2d(self.vor)
        fig.set_size_inches(12, 10)
        plt.xlabel('longitude',size=12)
        plt.ylabel('latitude',size=12)
        plt.title('Voronoi Diagram Colombia',size=12)
        plt.plot(Borders[:,1],Borders[:,0])
        plt.plot()
        
        
    def get_vertex(self,indices):
        #Función para tener las coordenadas de los vertices de cada region enumeradas 
        ver=[]
        vertices=self.vor.vertices
        for i in indices:
            if i!=[] and -1 not in indices:
                ver.append(vertices[i])
        return ver
    
    def Punto_2(self):
        #Tiempo de ejecución O(log n) Estoy calculando el area de cada diagramade voronoi 
        #el cual tiene un tiempo log n y lo estoy haciendo con todas las regiones
        
        Borders=self.borders.to_numpy()
        area=[]
        r=self.vor.regions
        point_region=self.vor.point_region
        n=len(point_region)
        for i in range(n):
            zona=r[point_region[i]]
            if -1 in zona or zona==[]:
                area.append(1)
            if -1 not in zona and zona!=[]:
                m=self.get_vertex(zona)
                area.append(ConvexHull(m).volume)
                
        print("Min index", area.index(min(area)), "Max index" ,area.index(max(area)))
        print("Area minima", self.airports['name'][area.index(min(area))], 'Area maxima', self.airports["name"][area.index(max(area))])
        
        minindex=point_region[area.index(min(area))]     
        maxindex=point_region[area.index(max(area)) ]
        
        
        indices=self.vor.regions #Indices de las areas a calcular
        ola=self.get_vertex(indices)
        
                
        fig=voronoi_plot_2d(self.vor)
        fig.set_size_inches(12, 10)
        plt.plot(Borders[:,1],Borders[:,0])
        plt.xlabel('longitude',size=12)
        plt.ylabel('latitude',size=12)
        plt.title('Max Area Yellow',size=12)
        # plt.plot(ola[20][:,0],ola[20][:,1], color='b')
        plt.fill(ola[maxindex-1][:,0],ola[maxindex-1][:,1], color='yellow')
        plt.plot()
        
        
        fig=voronoi_plot_2d(self.vor)
        fig.set_size_inches(12, 10)
        plt.plot(Borders[:,1],Borders[:,0])
        plt.xlabel('longitude',size=12)
        plt.ylabel('latitude',size=12)
        plt.title('Min Area Red',size=12)
        plt.plot(ola[minindex-1][:,0],ola[minindex-1][:,1], color='b')
        plt.fill(ola[minindex-1][:,0], ola[minindex-1][:,1], color='red')
        plt.plot()

    def Punto_3(self):
        #Tiempo de Ejecución O(n log n) Ya que estoy creando un arbol KdTree el cual tiene un tiempo de ejecución 
        #O(log n) y lo estoy usando n veces donde n es el len de la lista
        Borders=self.borders.to_numpy()
        l_1=[]
        for i in range(len(self.vor.vertices)):
            l_1.append(KDTree(self.vor.vertices).query_ball_point(self.vor.vertices[i], 1))
            
        k=min([(len(x)) for x in l_1])
        
        for i in l_1:
            if len(i)==k:
                p=i
                break
            
        circle1 = plt.Circle(self.vor.vertices[p[0]], 1, ec='r', fc='r')
        fig=voronoi_plot_2d(self.vor)
        fig.set_size_inches(12, 10)
        plt.xlabel('longitude',size=12)
        plt.ylabel('latitude',size=12)
        plt.title('New Airport Punto 3',size=12)
        plt.plot(Borders[:,1],Borders[:,0])
        plt.gcf().gca().add_artist(circle1)
        
        
    def Punto_4(self):
        #Calculamos el numero de vecinos para cada vertice y le agregamos 1 cada vez que encuentre un vecino
        
        vertices=self.vor.ridge_points
        vecinos=[0 for i in range(len(self.vor.points))]
        for i in vertices:
            vecinos[i[0]] += 1
            vecinos[i[1]] += 1
            
        min_index=vecinos.index(min(vecinos))
        max_index = vecinos.index(max(vecinos))
        name_min_airport=self.airports['name'][min_index]
        name_max_airport=self.airports['name'][max_index]
        
        print('Ciudad con más vecinos,' ,name_max_airport, 'Ciudad Con menos vecinos', name_min_airport)
        
        coord1=self.points[min_index]
        coord2=self.points[max_index]
        Borders=self.borders.to_numpy()
        #Graficamos las ciudades con mas y menos vecinos en un circulo para que pueda ser visible
                      
        circle2 = plt.Circle(coord1, 0.2, ec='r', fc='r')
        circle3=plt.Circle(coord2, 0.2, ec='r', fc='y')
        fig=voronoi_plot_2d(self.vor)
        fig.set_size_inches(12, 10)
        plt.xlabel('longitude',size=12)
        plt.ylabel('latitude',size=12)
        plt.title('Ciudades Con más veinos y menos(rojo, amarillo)',size=12)
        plt.plot(Borders[:,1],Borders[:,0])
        plt.gcf().gca().add_artist(circle2)
        plt.gcf().gca().add_artist(circle3)
        
                  
                

#Código adaptado de https://github.com/implse/Closest_Pair_Of_Points_Notebook/blob/master/Closest%20Points.ipynb
#Tiempo de ejecución O(n log n)

    
    def distance(self,p1, p2):
        d = math.sqrt(((p2[0] - p1[0])** 2) + ( (p2[1] - p1[1]) ** 2))
        return d
    
    ## Solo se utiliza cuando a n le quedan 3 elementos
    def closest_brute_force(self,points):
        min_dist = float("inf")
        p1 = None
        p2 = None
    
        for i in range(len(points)):
            for j in range(i+1, len(points)):
                d = self.distance(points[i], points[j])
    
                if d < min_dist:
                    p1 = points[i]
                    p2 = points[j]
                    min_dist = d
        return p1, p2, min_dist
    
    
    def recursive_closest(self,xsorted, ysorted):
        n = len(xsorted)
        if n <= 3:
            return self.closest_brute_force(xsorted)
        else:
            midpoint = xsorted[n//2]
            xsorted_left = xsorted[:n//2]
            xsorted_right = xsorted[n//2:]
            ysorted_left = []
            ysorted_right = []
            for point in ysorted:
                ysorted_left.append(point) if (point[0] <= midpoint[0]) else ysorted_right.append(point)
            (p1_left, p2_left, delta_left) = self.recursive_closest(xsorted_left, ysorted_left)
            (p1_right, p2_right, delta_right) = self.recursive_closest(xsorted_right, ysorted_right)
            (p1, p2, delta) = (p1_left, p2_left, delta_left) if (delta_left < delta_right) else (p1_right, p2_right, delta_right)
            in_band = [point for point in ysorted if midpoint[0]-delta < point[0] < midpoint[0]+delta]
            for i in range(len(in_band)):
                for j in range(i+1, min(i+7, len(in_band))):
                    d = self.distance(in_band[i], in_band[j])
                    if d < delta:

                        (p1, p2, delta) = (in_band[i], in_band[j], d)
            return p1, p2, delta
    
    
    def closest(self,points):
        xsorted = sorted(points, key=lambda point: point[0])
        ysorted = sorted(points, key=lambda point: point[1])
        return self.recursive_closest(xsorted, ysorted)
    
    # Function Call
    
    
    def Punto_5(self):
        
    
        closest_points = self.closest(self.points)
        Borders=self.borders.to_numpy()
        
        point_a, point_b, distance = closest_points
        print(point_a, point_b, 'Distancia entre los puntos ',distance)
        
        print('Aeropuertos a quitar, Tolemaida Air Base')
        print('Aeropuertos a quitar, Melgar Ab Airport')
        

        
        
        fig=voronoi_plot_2d(self.vor)
        fig.set_size_inches(12, 10)
        plt.plot(Borders[:,1],Borders[:,0])
        plt.plot(point_a,point_b, color='violet')
        plt.plot()
        
        
        x, y = zip(*self.points)
        x_a, y_a = point_a
        x_b, y_b = point_b
        
        plt.figure(figsize=(6, 6))
        plt.xlabel('longitude',size=12)
        plt.ylabel('latitude',size=12)
        plt.title('Closest Points',size=12)
        plt.scatter(x, y, color="grey")
        plt.scatter(x_a, y_a, color="red")
        plt.scatter(x_b, y_b, color="blue")
        plt.show()
        
        
        
        new_airport=[(point_a[0]+point_b[0])/2,(point_a[1]+point_b[1])/2]
        
        print('Nuevas cordenadas Aeropuerto', new_airport)
        plt.figure(figsize=(6, 6))
        plt.xlabel('longitude',size=12)
        plt.ylabel('latitude',size=12)
        plt.title('New Airport',size=12)
        plt.scatter(x, y, color="grey")
        plt.scatter(new_airport[0], new_airport[1], color="red")
        
        plt.show()
            