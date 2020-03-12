# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 17:59:43 2019

@author: Jorge A. Matias
"""
#se utiliza el algoritmo de la universidad politécnica de minas del documento de 
#Ultano Kindelan
#Módulo Compressed row storage
import numpy as np
def crs(matriz):
    ###CRS##
    #Al utilizar python se considera el indice inicial = 0
    x,y = matriz.shape#se obtienen las dimensiones de la matriz
    val = []#lista para almacenar los valores no nulos
    col_ind = []#lista para almacenar los índices de la columna de cada valor
    ren_elem = []#Lista para almacenar la informacion de los renglones de la matriz
    cont = 0#indice inicial del vector de renglones
    
    for i in range(x):
        
        ren_elem.append(cont)
        
        for j in range(y):
            
                if matriz[i,j] != 0:
                    
                    val.append(matriz[i,j])#se alnacenan los valores de la matriz
                    col_ind.append(j)#se almacenan los valores de la columna
                    cont = cont+1#incremento de los indices de renglones
   
    #se transforman las listas a arreglos de numpy
    valor = np.array(val)
    col = np.array(col_ind)
    ren = np.array(ren_elem)
    
    return(valor,col,ren)
    
def prodmat(valor,col,ren,vector):##Producto matricial utilizando CRS
    
    res=[]#lista para almacenar los resultados del producto
    
    for i in range(0,ren.size):
        suma = 0
        if i != ren.size-1:
            for j in range(ren[i],ren[i+1]):
                
                suma = suma + valor[j]*vector[col[j]]
            
            res.append(suma)
        else:
            for j in range(ren[i],valor.size):
                suma = suma + valor[j]*vector[col[j]]               
            res.append(suma)
    prod=np.array(res)
    return(prod) 
###############################################################################
###Producto matriz con vector
#def prodmat(matriz,vector):
#    res=[]
#    dmat=matriz.shape
#    for i in range(dmat[0]):
#        suma=0
#        for j in range(dmat[1]):
#            suma=suma+matriz[i,j]*vector[j]
#        res.append(suma)
#    prod=np.array(res)
#    return(prod)
###############################################################################
####Producto punto entre dos vectores
def dot(vec1,vec2):
    if vec1.shape == vec2.shape:#los vectores deben tener el mismo numero de elementos
        suma = 0
        for i in range(vec1.shape[0]):
            
            suma = suma + vec1[i]*vec2[i]
        res = np.array(suma)
    return(res)


def gradc(matriz,vector,tol):
    dim = matriz.shape[0]#tamaño de la matriz
    #Datos a CRS
    val,col_ind,ren_in = crs(matriz)#cambio a formato CRS
    #vectores iniciales
    x = np.ones(dim)#Solucion inicial, vector lleno de valor = 1
    r = np.copy(vector)-prodmat(val,col_ind,ren_in,x)#residuo
    d = np.copy(r)#direccion de descenso"d"
    beta = dot(r,r)
    maxitera = 1000
    i = 0
    while i<maxitera:
        print(i)
        z = prodmat(val,col_ind,ren_in,d)
        rho = beta/(dot(d,z))
        x = x+rho*d
        r = r-rho*z
        gamma = beta
        beta = dot(r,r)
        alfa = (beta/gamma)
        d = r + alfa*d
        #comprobacion
        if np.amax(r)<=tol:#norma del error maximo
            break
        i = i+1
        #print('iteracion',i)
    #print(x)   
    return(x)

tol = 0.000000001
mat = np.array([[4, -1, 0, -1, 0, 0, 0, 0, 0],
              [-1, 4, -1, 0, -1, 0, 0, 0, 0], 
              [0, -1, 4, 0, 0, -1, 0, 0, 0], 
              [-1, 0, 0, 4, -1, 0, -1, 0, 0], 
              [0, -1, 0, -1, 4, -1, 0, -1, 0], 
              [0, 0, -1, 0, -1, 4, 0, 0, -1],
              [0, 0, 0, -1, 0, 0, 4, -1, 0], 
              [0, 0, 0, 0, -1, 0, -1, 4, -1],
              [0, 0, 0, 0, 0, -1, 0, -1, 4]])
vec = np.array([25.0, 50.0, 150.0, 0, 0, 50.0, 0, 0, 25.0])

resultado = gradc(mat,vec,tol)