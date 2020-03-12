# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 15:27:55 2019

@author: Jorge Antonio Matías López
"""
#Metodos numéricos para resolver sistemas de ecuaciones nxn, Modulos Jacobi y gauss seidel
import numpy as np
import math
from time import time

'''
Supuestos:
Las matrices a resolver son simétricas positivas definidas, resultado de la discretización de 
ecuaciones diferenciales parciales. Todas reciben matrices dispersas, y en el caso de 
gauss-Seidel CRS, la matriz se convierte al formato CRS y despues se resuelven las operaciones.
Funciones: 
Resolvedores: Jacobi, Gauss-seidel, Gradiente conjugado
Operaciones matriciales: Productos matriz-Vector, producto punto


''' 
def jacobi(matriz,vect,tol):#Método de Jacobi
    print('Jacobi')
    dim = matriz.shape#tamaño de la matriz
    x0 = np.zeros(dim[0])#solucion inicial
    comp = np.zeros(dim[0])#arreglo vacío para la comprobación
    maxitera = 100#número máximo de iteraciones
    res = np.zeros(dim[0])#vector de resultados
    error = []#Lista vacía para evaluar el error
    k = 0
    t_ini = time()#evaluación de tiempo de computo(inicial)
    while k < maxitera:
        suma = 0
        k = k+1
        for ren in range(0, dim[0]):
            suma = 0
            for col in range(0, dim[1]):
                if (col != ren):
                    suma += matriz[ren,col] * x0[col]               
            res[ren] = (vect[ren] - suma)/matriz[ren, ren]
        del error[:]
            #Comprobación
        for ren in range(0, dim[0]):
            suma = 0
            for col in range(0, dim[1]):
                suma += matriz[ren,col]*res[col]  
                
            comp[ren] = suma
            
            dif=abs(comp[ren]-vect[ren])
            error.append(dif)#vector de errores
        
        for i in range(0, dim[0]):
            x0[i]=res[i]#actualización del resultado
        if all( i <= tol for i in error) == True:#condición para que termine si error<tolerancia
            break
        print('Iteracion',k)
    
    t_fin = time()
    t_ejecucion = round(t_fin-t_ini,20)
    print('El tiempo de ejecución Jacobi es '+str(t_ejecucion)+' segundos')
    return(res)
###############################################################################    
def gauss(matriz,vect,tol):#Método de Gauss-Seidel
    
    print('Gauss-Seidel')
    
    dim = matriz.shape
    comp = np.zeros(dim[0])
    maxitera = 100
    res = np.zeros(dim[0])
    error = []
    k = 0
    t_ini = time()
    while k < maxitera:
        suma = 0
        k += 1
        for ren in range(0, dim[0]):
            suma = 0
            for col in range(0, dim[1]):
                if (col != ren):
                    suma += matriz[ren,col]*res[col]               
            res[ren] = (vect[ren] - suma)/matriz[ren, ren]
        del error[:]
            #Comprobación
        for ren in range(0,dim[0]):
            suma = 0
            for col in range(0,dim[1]):
                suma += matriz[ren,col]*res[col]  
            comp[ren] = suma
            dif = abs(comp[ren]-vect[ren])
            
            error.append(dif)
#            print('Error en x"',ren,'"=',error[ren])
        if all( i<= tol for i in error) == True:
            break
        print('Iteracion',k)
    t_fin = time()
    t_ejecucion = round(t_fin-t_ini,20)
    print('El tiempo de ejecución gauss es '+str(t_ejecucion)+' segundos')
    return(res)
###############################################################################
def gauss_crs(matriz,vector,tol):#Gauss-Seidel usando CRS
    ###CRS
    x,y = matriz.shape
    val = []#lista para almacenar los valores no nulos
    col_ind = []#Lista para almacenar el índice de columna de cada dato
    ren_elem = []#Lista para almacenar los datos de los renglones
    di = []#lista para almacenar la diagonal
    cont = 0
    for i in range(x):
        ren_elem.append(cont)
        for j in range(y):
            if i != j:
                if matriz[i,j] != 0:
                    val.append(matriz[i,j])
                    col_ind.append(j)
                    cont += 1
            else:
                di.append(matriz[i,j])
    #Transforman las listas en arreglos de numpy           
    valor = np.array(val)
    col = np.array(col_ind)
    ren = np.array(ren_elem)
    diag = np.array(di)
   ##GaussSeidel 
    maxitera = 100
    res = np.zeros(x)#arreglo inicial de resultados
    exa = np.linalg.solve(matriz,vector)
    error = []#lista para evaluar el error
    k = 0
    t_ini = time()
    while k < maxitera:
        suma = 0
        k += 1
        for i in range(0,ren.size):
            suma = 0
            if i != ren.size-1:
                for j in range(ren[i],ren[i+1]):
                    suma += valor[j] * res[col[j]]  
                    
                res[i] = (vector[i] - suma)/diag[i]
            
            else:
                for j in range(ren[i],valor.size):
                    suma += valor[j]*res[col[j]]               
                res[i] = (vector[i]-suma)/diag[i]
    
        del error[:]#Elimina los valores de la lista en cada iteracion
        
        
            #Comprobación
        for i in range(0,res.size):
            dif = abs(exa[i]-res[i])
            error.append(dif)#vector de errores
        if all( i <= tol for i in error) == True:
            break
        print('Iteracion',k)
    t_fin = time()
    t_ejecucion = round(t_fin-t_ini,10)
    print('El tiempo de ejecución CRS_gauss es '+str(t_ejecucion)+' segundos')
    return(res)

###Gradiente conjugado sin CRS#################################################
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
    if vec1.shape == vec2.shape:
        suma = 0
        for i in range(vec1.shape[0]):
            suma += vec1[i]*vec2[i]
        res=np.array(suma)
    return(res)

def norma(r):#Norma del vector error
    n = r.size
    suma = 0
    for i in range(n):
        #error
        suma += r[i]**2
    error = math.sqrt(suma)
    return(error)
#Metodo de Gradiente conjugado
def gradc(matriz,vector,tol):
    dim = matriz.shape[0]#tamaño de la matriz
    #vector inicial
    x = np.ones(dim)
    r = np.copy(vector)-prodmat(matriz,x)#residuo
    d = np.copy(r)#direccion de descenso"d"
    beta = dot(r,r)
    maxitera = 100
    i = 0
    t_ini = time()
    while i < maxitera:
        z = prodmat(matriz,d)
        rho = beta/(dot(d,z))
        x = x+rho*d
        r = r-rho*z
        gamma = beta
        beta = dot(r,r)
        alfa = (beta/gamma)
        d = r + alfa*d
        if norma(r) <= tol:
            break
        i += 1
        #print('iteracion',i)
    t_fin = time()
    t_ejecucion = round(t_fin-t_ini,10)
    #print(x)
    print('El tiempo de ejecución Gradiente conjugado es '+str(t_ejecucion)+' segundos')    
    return(x)



mat=np.array([[10,-1,2,0],
              [-1,11,-1,3],
              [2,-1,10,-1],
              [0,3,-1,8]],dtype=np.float64)
vector=np.array([6,25,-11,15],dtype=np.float64)
toler=0.0000001
soljacobi=jacobi(mat,vector,toler)
solgauss=gauss(mat,vector,toler)
sollinal=np.linalg.solve(mat,vector)
solgauss_crs=gauss_crs(mat,vector,toler)
sol_gradc=gradc(mat,vector,toler)