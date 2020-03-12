# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 17:58:15 2019

@author: jor_a
"""
#Modulo gauss-jordan
#matriz nxn
#vector tamaño n
import numpy as np

def gauss_j(matr,vect):
    matriz=matr.copy()
    print(matriz is mat)
    vector=vect.copy()
    print(vector is vec)
    x,y=matriz.shape
    for i in range(x):
        j=i
        if matriz[i,j]!=1:
            vector[i]=vector[i]/matriz[i,j]
            matriz[i,:]=matriz[i,:]/matriz[i,j]
        for ren in range(i+1,x):
            if matriz[ren,j]!=0:
                a=-matriz[ren,j]
                vector[ren]=vector[i]*a+vector[ren]
                for col in range(j,x):
                    matriz[ren,col]=matriz[i,col]*a+matriz[ren,col]
            
    
    return(matriz,vector)

mat=np.array([[2,1,3],
              [3,-1,-2],
              [0,-1,1]],dtype=np.float64)
    
vec=np.array([7,5,10],dtype=np.float64)

a,b=gauss_j(mat,vec)
sol=np.linalg.solve(mat,vec)