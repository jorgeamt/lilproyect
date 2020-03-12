# -*- coding: utf-8 -*-
"""
Created on Mon Jun 10 14:21:23 2019

@author: jor_a
"""
#GMRES
#Datos de entrada=Matriz A(nxn), constante b, 
'''
Hay que optimizar el algoritmo gmres, se puede intentar usar minimos cuadrados en el sistema interno

'''



import numpy as np
import math
from time import time
###############################################################################    
def gauss_j(matr,vect):
    t_inig = time()
    matriz = matr.copy()
    vector = vect.copy()
    x,y = matriz.shape
    for i in range(x):
        j = i
        if matriz[i,j] != 1:
            vector[i] = vector[i]/matriz[i,j]
            matriz[i,:] = matriz[i,:]/matriz[i,j]
        for ren in range(i+1,x):
            if matriz[ren,j] != 0:
                a = -matriz[ren,j]
                vector[ren] = vector[i]*a + vector[ren]
                for col in range(j,x):
                    matriz[ren,col] = matriz[i,col]*a + matriz[ren,col]
    for j in range(x-1,-1,-1):
        i = j
        for ren in range(j-1,-1,-1):
            if matriz[ren,j] != 0:
                a = -matriz[ren,j]
                vector[ren] = vector[i]*a + vector[ren]
                for col in range(j,j-1,-1):
                    matriz[ren,col] = matriz[i,col]*a + matriz[ren,col]
    t_fing = time()
    t_ejecucion = round(t_fing - t_inig, 20)
    print('El tiempo de ejecución gauss_j es '+str(t_ejecucion)+' segundos')
    return(vector)
    
###############################################################################
def norma(r):#Norma del vector 
    n = r.size
    suma = 0
    for i in range(n):
            #error
            suma = suma + r[i]**2
    error = math.sqrt(suma)
    return(error)
def prodmat(matriz,vector):
    res=[]
    dmat = matriz.shape
    for i in range(dmat[0]):
        suma = 0
        for j in range(dmat[1]):
            suma = suma + matriz[i,j]*vector[j]
        res.append(suma)
    prod = np.array(res)
    return(prod)
###############################################################################
####Producto punto entre dos vectores
def dot(vec1,vec2):
    if vec1.shape == vec2.shape:
        suma = 0
        for i in range(vec1.shape[0]):
            suma = suma + vec1[i]*vec2[i]
        res = np.array(suma)
    return(res)
def gmres(mat,vec):
    t_ini = time()
    t_fin=time()
    x,y = mat.shape
    x0 = np.zeros(x)
    m = x
    h = np.zeros([m+1,m])
    V = np.zeros([m,m+1])
    r = vec-prodmat(mat,x0)
    normr = norma(r)
    v = r/normr
    V[:,0] = v
    for j in range(m):
        vj1 = np.zeros(x)
        vj = V[:,j]
        Avj = prodmat(mat,vj)
        for i in range(j+1):
            vi = V[:,i]
            h[i,j] = dot(Avj,vi)
        for i in range(j+1):
            vi = V[:,i]
            vj1 = vj1 + h[i,j]*vi
        vj1 = -vj1 + Avj
        #Se obtiene V y H 
        h[j+1,j] = norma(vj1)
        V[:,j+1] = vj1/h[j+1,j]
    h = h[0:m,:]
    beta = np.zeros(m)
    beta[0] = normr
    print(beta)


    ykl = np.linalg.solve(h,beta)
    print(ykl)
    yk = gauss_j(h,beta)
    print(yk)
    for i in range(m):
        x0 = x0 + yk[i]*V[:,i]
    t_fin=time()
    t_ejecucion=round(t_fin-t_ini,20)
    print('El tiempo de ejecución gmres es '+str(t_ejecucion)+' segundos')
    return (V,h,beta,x0)
'''mat = np.array([[2,1,-1,3],
              [1,-2,2,1],
              [3,-2,1,4],
              [2,1,2,8]],dtype=np.float64)

vector = np.array([8,6,2,4],dtype=np.float64)
'''
mat = np.random.rand(300,300)
vector = np.random.rand(300)
V,h,beta,x = gmres(mat,vector)

sol_gauss=gauss_j(mat,vector)

sol = np.linalg.solve(mat,vector)
#a=np.linalg.inv(A)



