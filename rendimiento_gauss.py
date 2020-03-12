# -*- coding: utf-8 -*-
"""
Created on Tue Jul 2 17:01:35 2019

@author: Jorge Antonio Matías López
"""

#prueba gauss jordan
import numpy as np
import matplotlib.pyplot as plt
from time import time


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
    return(vector,t_ejecucion)
"""
tiempo = []
incognitas = []
for k in range(6,606,6):
    mat = np.random.rand(k,k)
    vector = np.random.rand(k)

    sol_gauss,t=gauss_j(mat,vector)   
    tiempo.append(t)
    incognitas.append(k)
Y1=np.array(tiempo)
X1=np.array(incognitas)

f, (ax1) = plt.subplots()
ax1.plot(X1,Y1)
ax1.grid()
ax1.set(xlabel = 'Tiempo de ejecución', ylabel = 'Número de incógnitas',
       title='Desempeño') 
plt.show()
"""
mat = np.random.rand(50,50)
vector = np.random.rand(50)
sol_gauss,t=gauss_j(mat,vector)