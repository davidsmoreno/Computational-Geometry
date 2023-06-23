#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  2 08:53:37 2021

@author: david
"""

#Crear clase y constructor de puntos con su elevación
#Hacer una triangulación de esos puntos

import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import Delaunay
import matplotlib.tri as mtri
from scipy.spatial import Delaunay
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import minimum_spanning_tree
import networkx as nx

from scipy.interpolate import LinearNDInterpolator
import pandas as pd







class TIN:
    #atributos
    #num_points: Número de puntos
    #triangulation: scypy.spatial.Delunay
    #elevations: list of floats(one-to one )
    
    
    def __init__(self,num_points,puntos,z):
        self.num_points=num_points
        self.puntos=puntos
        self.z=z
        self.triangulacion=Delaunay(self.puntos)
    
    def graph_punto1_2d(self):
        plt.tripcolor(self.puntos[:,0], self.puntos[:,1],self.triangulacion.simplices,self.z,cmap='coolwarm', edgecolors='k')
        plt.show()
        
    def graph_punto1_3d(self):
        triang = mtri.Triangulation(self.puntos[:,0], self.puntos[:,1], triangles=self.triangulacion.simplices)
        fig, ax = plt.subplots(subplot_kw =dict(projection="3d"))
        ax.plot_trisurf(triang, self.z, cmap = 'coolwarm')
        plt.get_cmap('coolwarm')

        plt.show()
        
    def Interpolacion_punto2(self,puntoxy):
        
        t_indice=self.triangulacion.find_simplex(puntoxy) #Calcula en qu'e triangulo est'a 
        triangulo=self.triangulacion.simplices[t_indice] #Nos da el indice de los puntos que conectan el triangulo
        elev=self.triangulacion.points[triangulo] # Puntos en 2d que conectan el triangulo
        indices_z=[i for i in triangulo] # Altura de los puntos en 2d
    
    
        f=LinearNDInterpolator(elev,self.z[indices_z]) #hacemos la interporlacion lineal con los 3 puntos en 3d del triangulo

        print(f(puntoxy)) #imprimimos la altura del punto a partir de la interporlacion lineal
    def gratafica_Interpolacion_punto2(self,puntoxy):
        t_indice=self.triangulacion.find_simplex(puntoxy) #Calcula en qu'e triangulo est'a 
        triangulo=self.triangulacion.simplices[t_indice] #Nos da el indice de los puntos que conectan el triangulo
        elev=self.triangulacion.points[triangulo] # Puntos en 2d que conectan el triangulo
        plt.figure()
        plt.scatter(elev[:,0],elev[:,1])
        plt.scatter(puntoxy[0],puntoxy[1])
        plt.plot(elev[:,0],elev[:,1])
        plt.show()
        
    def Area_Cuadrilatero(self,puntos):

        return 0.5 * np.abs(np.dot(puntos[:,0], np.roll(puntos[:,1], 1)) - np.dot(puntos[:,1], np.roll(puntos[:,0], 1)))
    
    def Largest_area(self,puntoxy):
        
        Area=[]
        t_indice=self.triangulacion.find_simplex(puntoxy) #Calcula en qu'e triangulo est'a 
        triangulo=self.triangulacion.simplices[t_indice]
        punt_triangulo=self.triangulacion.points[triangulo] 
        Vecinos_indices=self.triangulacion.neighbors

        vecinos_triangulo=Vecinos_indices[t_indice]
    
    
        punt_vecinos=self.triangulacion.points[vecinos_triangulo]
    
        for i in range(len(punt_vecinos)):
            Area.append(self.Area_Cuadrilatero(np.concatenate((punt_triangulo, [punt_vecinos[i]]),axis=0)))
    
        punt_areamax=max(range(len(Area)), key=Area.__getitem__)
        puntos=np.concatenate((punt_triangulo,[punt_vecinos[punt_areamax]]),axis=0)
        return puntos
    
    
    def grafica_Largest_area(self,puntos):
        plt.scatter(puntos[:,0],puntos[:,1])
        plt.plot(puntos[:,0], puntos[:,1])
        plt.show()
        
        
    def Angulo_Triangulo(self,puntos):
        #Adaptación https://stackoverflow.com/questions/5122372/angle-between-points
        puntos=np.array(puntos)
        angulo=[]
        a=puntos[2]-puntos[0]
        b=puntos[1]-puntos[0]
        c=puntos[2]-puntos[1]
        
        for i,j in ((a,b),(a,c),(b,-c)):
            p=np.dot(i,j)
            n=np.linalg.norm(i) * np.linalg.norm(j)
            angulo.append(np.arccos(p/n)*180/np.pi)
        return angulo
        
        
    
    
    def angulos_tringulacion(self):
        angulos=[]
        for i in range(len(self.triangulacion.simplices)):
            triangulo=self.triangulacion.simplices[i]
            punt_triangulo=self.triangulacion.points[triangulo] 
            angulo=self.Angulo_Triangulo(punt_triangulo)
            angulos.append(angulo)
            
            B = np.reshape(angulos, (-1, 3))
        return B 
    
    
    def max_min_angulos(self):
        B=self.angulos_tringulacion()
        maximo = max([max(l) for l in B])
        minimo = min([min(l) for l in B])
        #Triangulo mayor y menor angulo
        print("Ratio of the maximum angle and the minimum", maximo/minimo)
        
        
    def SpanningTree(self):
        
        z1=self.z.reshape(self.num_points,1)
        U=np.concatenate((self.puntos,z1),axis=1)

        #Create a Adyacent Weigth Matrix fuente https://www.tutorialfor.com/questions-308657.htm
        df=pd.DataFrame(U,columns=["node_1","node_2","weight"],)
        name_to_node={name: i for i, name in enumerate(np.unique(df[["node_1", "node_2"]].values))}
        n_nodes = len(name_to_node)
        A = np.zeros((n_nodes, n_nodes))
        for row in df.itertuples():
            n1 = name_to_node[row.node_1]
            n2 = name_to_node[row.node_2]
            A[n1, n2] += row.weight
            A[n2, n1] += row.weight
            
        Matrizadyacencia=A
        X = csr_matrix(Matrizadyacencia)
        Tcsr = minimum_spanning_tree(X)
        Minimu_tree=Tcsr.toarray().astype(int)
        return Minimu_tree
            
           

          
        
        
        


        
