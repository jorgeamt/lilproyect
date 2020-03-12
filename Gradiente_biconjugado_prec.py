# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 11:20:47 2019

@author: Jorge Antonio Matías López
"""

#Gradiente biconjugado precondicionado

import numpy as np
from time import time
#import matplotlib.pyplot as plt

def maxerror(r):
    val = abs(r[0])
    for i in r:
        if abs(i) >= val:
            val = abs(i)

    return(val)

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
                
                suma += valor[j]*vector[col[j]]
            
            res.append(suma)
        else:
            for j in range(ren[i],valor.size):
                suma += valor[j]*vector[col[j]]               
            res.append(suma)
    prod = np.array(res)
    
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
            
            suma += vec1[i]*vec2[i]
        res = np.array(suma)
    return(res)


def gradbic(matriz,vector,tol):
    t_inig = time()
    
    dim = matriz.shape[0]#tamaño de la matriz
    
    #Datos a CRS
    val,col_ind,ren_in = crs(matriz)#cambio a formato CRS
    At = np.transpose(matriz)#Transpuesta de la matriz
    valt,col_indt,ren_int = crs(At)#Matriz transpuesta convertida a CRS
    
    #precondicionador
    #se usa la pseudo inversa, es decir la inversa de la diagonal
    #se tiene entonces que M = transpuesta de M
    # En caso de usar otro precondicionador se deberá calcular la transpuesta y su formato crs  
    M = np.identity(dim)  
    for i in range(dim):
        for j in range(dim):
            if i == j:
                M[i,j] *= 1/matriz[i,j]
    Mt = np.transpose(M)
    
    valm, col_indm, ren_inm = crs(M)
    
    valmt, col_indmt, ren_inmt = crs(Mt)
    
    #vectores iniciales 
    x = np.ones(dim)#Solucion inicial, vector lleno de valor = 1     
    r = np.copy(vector)-prodmat(val,col_ind,ren_in,x)#residuo(Gradiente inicial)
    r_ps = np.copy(r)# Pseudo gradiente inicial
    
    q = prodmat(valm, col_indm, ren_inm,r)
    q_ps = prodmat(valmt, col_indmt, ren_inmt, r_ps)
    
    p = np.copy(q)#dirección de descenso inicial
    p_ps = np.copy(q_ps)
    
    maxitera = 1000
    i = 0
   
    while i < maxitera:
        
        beta = dot(q,r_ps)
        
        w = prodmat(val,col_ind,ren_in,p)
        
        w_ps = prodmat(valt,col_indt,ren_int,p_ps)
        
        
        rho = beta/(dot(w, p_ps))
        
        x += rho * p
        r -= rho * w
        r_ps -= rho * w_ps
        
        q = prodmat(valm, col_indm, ren_inm,r)        
        q_ps = prodmat(valmt, col_indmt, ren_inmt, r_ps)
                
        gamma = beta
        beta = dot(r_ps,q)
        alfa = (beta/gamma)
        
        p = q + alfa * p
        p_ps = q_ps + alfa * p_ps
        
        #comprobacion
        
        if maxerror(r)<=tol:#norma del error maximo
            
            print(tol)
            print(maxerror(r))
            break
        i += 1
        
    t_fing = time()
    print('Iteraciónes:',i)
    t_ejecucion = round(t_fing - t_inig, 20)
    print('El tiempo de ejecución Gradiente biconjugado es '+str(t_ejecucion)+' segundos')
    return(x,t_ejecucion)
    
tol = 0.000001
mat = np.random.rand(50,50)
vector = np.random.rand(50)
sol_numpy = np.linalg.solve(mat,vector)
sol_gradb, t2 = gradbic(mat,vector,tol) 
