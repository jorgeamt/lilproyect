# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 09:33:05 2019

@author: jor_a
"""

#matriz nxn
#vector tamaño n
import numpy as np
def gauss_j(matr,vect):
    '''
    En esta función utiliza dos argumentos, el primero es la matriz cuadrada de tamaño nxn
    y el segundo argumento es el vector de datos independientes de tamaño n
    '''
    #Se crean copias de la matriz y el vector para no modificar los originales
    matriz = matr.copy()
    vector = vect.copy()
    x,y = matriz.shape # Obtenemos las dimensiones de la matriz
    
    
    for i in range(x):
        # como primer paso verificamos que el valor de la diagonal en cada renglón sea 1
        j = i
        
        if matriz[i,j] != 1:
            
            vector[i] = vector[i]/matriz[i,j]
            
            matriz[i,:] = matriz[i,:]/matriz[i,j]

        #los elementos de la columna que se encuentran debajo de la diagonal deben ser cero
        
        for ren in range(i + 1, x):#Recorremos los renglones hacia abajo
            
            if matriz[ren,j] != 0:
                
                a = -matriz[ren,j]# Guardamos el valor diferente de cero con signo contrario
                #se hace la multiplicación del valor guardado y se suma con los renglones inferiores en el vector
                vector[ren] = vector[i] * a + vector[ren] 
                
                #se hace la multiplicación del valor guardado y se suma con los renglones inferiores en la matriz
                for col in range(j, x):
                    
                    matriz[ren,col] = matriz[i, col]*a + matriz[ren,col] 
            
    #En este punto ya tenemos la matriz triangular superior, para terminar el proceso de Gauss-Jordan 
    #debemos realizar el último proceso de manera inversa
    for j in range(x-1, -1, -1):
        i = j
        for ren in range(j-1, -1, -1):
            if matriz[ren,j] != 0:
                
                a = -matriz[ren,j]
                vector[ren] = vector[i]*a + vector[ren]
                
                for col in range(j,j-1,-1):
                    matriz[ren,col]=matriz[i,col]*a+matriz[ren,col]

    return(matriz,vector)#devolvemos la matriz identidad y la solución del sistema de ecuaciones

mat=np.array([[2,1,3],
              [3,-1,-2],
              [0,-1,1]],dtype=np.float64)
    
vec=np.array([7,5,10],dtype=np.float64)
a,b = gauss_j(mat,vec)

print(a)
print(b)