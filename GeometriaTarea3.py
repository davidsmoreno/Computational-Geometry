#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 27 01:08:39 2021

@author: david
"""

import numpy as np
import random 
from scipy import spatial
from scipy.spatial import KDTree
import matplotlib.pyplot as plt 
from matplotlib.patches import Rectangle


class FlightOperations:
    
        
    def __init__(self,points_x,points_y):
        
        self.data=np.dstack([points_x.ravel(),points_y.ravel()])[0]
        self.n_points=len(self.data)
        self.Arbol=spatial.KDTree(self.data,leafsize=1)
        self.points_x=points_x
        self.points_y=points_y
        


    def Punto_1(self,point):
        distancia, indice = self.Arbol.query(point)
        barco=self.Arbol.data[indice] ## Barco mas cercano 
        fig = plt.figure()
        fig.set_size_inches(5, 5)
        plt.xlabel('Eje x',size=10)
        plt.ylabel('Eje y',size=10 )
        plt.title('Punto 1',size=10)
        plt.scatter(self.data[:,0],self.data[:,1],marker = '>',label='Barcos')
        plt.scatter(point[0],point[1],marker = 'o', color='g',label='vessel')
        plt.scatter(barco[0],barco[1],marker = '>',color='r',label='Barco más cercano')
        plt.xlim(-1,1)
        plt.ylim(-1,1)
        plt.legend(bbox_to_anchor=(1.5, 1.05))
        plt.plot()
        
        print('La Coordenada Vessel(Punto) es ',point, 'La coordenada del barco más cercano es ', barco)
        
    def Punto_2(self,r):
        barco1=self.Arbol.data[0]
        indice2 =self.Arbol.query_ball_point(barco1,r=r/2,p=np.inf)
        barco2=self.Arbol.data[indice2] ## Barco mas cercano 
        fig = plt.figure()
        fig.set_size_inches(5, 5)
        plt.xlabel('Eje x',size=10)
        plt.ylabel('Eje y',size=10)
        plt.title('Punto dos ',size=10)
        plt.scatter(self.data[:,0],self.data[:,1],marker = '>',label='Barcos')
        plt.scatter(barco2[:,0],barco2[:,1],marker='>',color='r',label='Barcos en el rango')
        plt.scatter(barco1[0],barco1[1],marker='>',color='k',label='Barco Elegido')
        plt.gca().add_patch(Rectangle(barco1-r/2,r,r,linewidth=1,edgecolor='r',facecolor='none'))
        plt.xlim(-1,1)
        plt.ylim(-1,1)
        plt.legend(bbox_to_anchor=(1.5, 1.05))
        plt.plot()
        
        print('La Coordenada del barco elegido es ',barco1)
        print('Coordenadas de los barcos Que están en el cuadrado de radio',r, 'Son :  ',barco2)


    def punto_3(self):
        Punto_oeste_index=np.argwhere(self.points_x==self.Arbol.maxes[0]).reshape(-1)
        Punto_norte_index=np.argwhere(self.points_y==self.Arbol.maxes[1]).reshape(-1)
        Punto_sur_index=np.argwhere(self.points_y==self.Arbol.mins[1]).reshape(-1)
        Punto_este_index=np.argwhere(self.points_x==self.Arbol.mins[0]).reshape(-1)
        
        b_oeste=self.Arbol.data[Punto_oeste_index]
        b_norte=self.Arbol.data[Punto_norte_index]
        b_sur=self.Arbol.data[Punto_sur_index]
        b_este=self.Arbol.data[Punto_este_index]
        
        fig = plt.figure()
        fig.set_size_inches(5, 5)
        plt.xlabel('Eje x',size=10)
        plt.ylabel('Eje y',size=10)
        plt.title('Punto 3',size=10)
        plt.scatter(b_oeste[0][0],b_oeste[0][1],marker = '>',color='y',label='Barco Oeste')
        plt.scatter(b_norte[0][0],b_norte[0][1],marker = '>',color='r',label='Barco Norte')
        plt.scatter(b_sur[0][0],b_sur[0][1],marker = '>',color='k',label='Barco sur')
        plt.scatter(b_este[0][0],b_este[0][1],marker = '>',color='g',label='Barco Este')
        plt.xlim(-1-0.1,1+0.1)
        plt.ylim(-1-0.1,1+0.1)
        plt.legend(bbox_to_anchor=(1.5, 1.05))
        plt.plot()
        
        print('Coordenadas Barco Oeste',b_oeste)
        print('Coordenadas Barco Norte',b_norte)
        print('Coordenadas Barco Sur',b_sur)
        print('Coordenadas Barco Este',b_este, 'Si dos puntos tiene coordenadas iguales queda el color verde, negro, rojo y amarillo en ese orden por encima de los demás')
        
        
 
                   
##Codigo adaptado de https://www.astroml.org/book_figures/chapter2/fig_kdtree_example.html

    def Plotkdtree(self,Nodo,minimos,maximos):
        if type(Nodo) is KDTree.leafnode: #Caso base, cuando el nodo sea la hoja, hemos terminado
            return None
        else:
            Min=minimos.copy()
            Max=maximos.copy()
            Max[Nodo.split_dim] = Nodo.split ##En la documentación vemos que existe una subclase, node donde
                                             ## Tenemos la propiedad de que podemos dividir un nodo en dos
            Min_r = minimos.copy()
            Min_r[Nodo.split_dim] = Nodo.split 
            Max_r= maximos.copy()   
            
            if Nodo.split_dim==0:  ##Si la dimención del nodo es 0, esto es. Estamos en la dirección x
                direccion_x=[Nodo.split,Nodo.split]  
                direccion_y=[minimos[1],maximos[1]]  ##Nos movemos en dirección paralela al eje x y buscamos el Punto y con menor y mayor coordenada
                color='r'                            ## Le agregamos el color azul
            else:
                direccion_x=[minimos[0],maximos[0]]
                direccion_y=[Nodo.split,Nodo.split]
                color='b'
            
    
            plt.plot(direccion_x,direccion_y,color) ##Vamos graficando las lineas con respecto a los puntos de mayor y menor coordenada
            plt.xlim(-1,1)
            plt.ylim(-1,1)
                                
            self.Plotkdtree(Nodo.less, Min, Max) ##LLamamos a la función recursiva Quitando el nodo más pequeno que dibujamos, hasta que lleguemos a 0 nodos
            self.Plotkdtree(Nodo.greater, Min_r, Max_r)  
        
            
    def Plot_Kdtree(self):
        fig = plt.figure()
        fig.set_size_inches(5, 5)
        plt.xlabel('Eje x',size=10)
        plt.ylabel('Eje y',size=10)
        plt.title('Punto 4',size=10)
        Nodo=self.Arbol.tree
        plt.scatter(self.Arbol.data[:,0],self.Arbol.data[:,1],marker = '>',label=None)
        minimos = np.copy(self.Arbol.mins)
        maximos = np.copy(self.Arbol.maxes)
        self.Plotkdtree(Nodo, minimos, maximos)
        plt.xlim(-1,1)
        plt.ylim(-1,1)
        plt.show()
        