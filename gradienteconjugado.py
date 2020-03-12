# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 13:56:09 2019

@author: jor_a
"""
#Gradiente conjugado
import numpy as np
import math
###Gradiente conjugado sin CRS
###############################################################################
###Producto matriz con vector
def prodmat(matriz,vector):
    res=[]
    dmat=matriz.shape
    for i in range(dmat[0]):
        suma=0
        for j in range(dmat[1]):
            suma=suma+matriz[i,j]*vector[j]
        res.append(suma)
    prod=np.array(res)
    return(prod)
###############################################################################
####Producto punto entre dos vectores
def dot(vec1,vec2):
    if vec1.shape==vec2.shape:
        suma=0
        for i in range(vec1.shape[0]):
            suma=suma+vec1[i]*vec2[i]
        res=np.array(suma)
    return(res)

def error(r):#Norma del vector error
    n=r.size
    suma=0
    for i in range(n):
            #error
            suma=suma+r[i]**2
    error=math.sqrt(suma)
    return(error)
def gradc(matriz,vector,tol):
    dim=matriz.shape[0]#tamaño de la matriz
    #vectores iniciales
    x=np.ones(dim)
    r=np.copy(vec)-prodmat(matriz,x)#residuo
    d=np.copy(r)#direccion de descenso"d"
    beta=dot(r,r)
    maxitera=10
    i=0
    while i<maxitera:
        z=prodmat(matriz,d)
        rho=beta/(dot(d,z))
        x=x+rho*d
        r=r-rho*z
        gamma=beta
        beta=dot(r,r)
        alfa=(beta/gamma)
        d=r+alfa*d
        print('x=',x,'rho=',rho,'r=',r,'alfa=',alfa,'d=',d)
        if error(r)<=tol:
            break
        i=i+1
        print('iteracion',i)
    return(x)
tol=0.000001
#mat=np.array([[4,3,0],
#              [3,4,-1],
#              [0,-1,4]])
#vec=np.array([24,30,-24])
mat=np.array([[4, -1, 0, -1, 0, 0, 0, 0, 0],
              [-1, 4, -1, 0, -1, 0, 0, 0, 0], 
              [0, -1, 4, 0, 0, -1, 0, 0, 0], 
              [-1, 0, 0, 4, -1, 0, -1, 0, 0], 
              [0, -1, 0, -1, 4, -1, 0, -1, 0], 
              [0, 0, -1, 0, -1, 4, 0, 0, -1],
              [0, 0, 0, -1, 0, 0, 4, -1, 0], 
              [0, 0, 0, 0, -1, 0, -1, 4, -1],
              [0, 0, 0, 0, 0, -1, 0, -1, 4]])
vec=np.array([25.0, 50.0, 150.0, 0, 0, 50.0, 0, 0, 25.0])
solucion=gradc(mat,vec,tol)
solucion=gradc(mat,vec,tol)
