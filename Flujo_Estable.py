# -*- coding: utf-8 -*-
"""
Created on Dom Jan 20 12:33:55 2019

@author: Jorge Antonio Matias López
"""
#Ejemplo de la solucion de un flujo estacionario mediante diferencias finitas con las condiciones:
#dx&dy=1
#S1=100, S2=0
#h0=50
#malla=[100, 100, 100, 100, 100],[50,50,50,50,50],[50,50,50,50,50],[50,50,50,50,50],[0,50,50,50,50]
#Ejemplo sencillo para utilizar graficas 
import numpy as np
import matplotlib.pyplot as plt
import math
def pltcont(matriz):#Función para invertir el orden de las filas de la matriz
    arreglo=[]
    for i in range(len(matriz)):
        arreglo.append(matriz[len(matriz)-1-i])
    return(arreglo)

malla=np.array([[100, 100, 100, 100],[50,50,50,50],[50,50,50,50],[0,50,50,50]], dtype=np.float64)
rencol=malla.shape
err=np.zeros(rencol)#matriz de error relativo
solc=np.zeros(rencol)#solucion conocida
x=[]
y=[]

ciclos=34#iteraciones

cont=1
while cont<=ciclos:
    
    for i in range(malla.shape[0]):
        
        for j in range(malla.shape[1]):
            
            solc[i][j]=malla[i][j]
            
            if len(malla)-1>i>0:#condición de frontera S1
                if j==0:#C. F. izquierda
                    malla[i][j]=round(((2*malla[i][j+1]++malla[i+1][j]+malla[i-1][j])/4),2)                    
                elif j==len(malla[i])-1:#C.F. derecha
                    malla[i][j]=round(((2*malla[i][j-1]++malla[i+1][j]+malla[i-1][j])/4),2)                    
                else:
                    malla[i][j]=round(((malla[i][j+1]+malla[i][j-1]+malla[i+1][j]+malla[i-1][j])/4),2)   
                
            elif i==len(malla)-1:#frontera inferior
                if j==0:#C. F. izquierda
                    pass#condicion de frontera S2                      
                elif j==len(malla[i])-1:#C.F. derecha
                    malla[i][j]=round(((2*malla[i][j-1]+2*malla[i-1][j])/4),2)
                    cont=cont+1                    
                else:
                    malla[i][j]=round(((malla[i][j+1]+malla[i][j-1]+2*malla[i-1][j])/4),2)
            err[i][j]=math.fabs(solc[i][j]-malla[i][j])#error absoluto            


X=np.array(x)
Y=np.array(y)

#graficas
calor=np.array(malla)#arreglo de numpy para graficar
contour=np.array(pltcont(malla))#arreglo de numpy ordenado

fig, ax = plt.subplots()#Solo genera una grafica
#im = ax.imshow(calor)
cs = ax.contourf(contour)#función para generar grafico
for i in range(len(malla)):
    for j in range(len(malla[i])):
        text = ax.text(j, i, contour[i, j],ha="center",
                       va="center", color="w")
ax.set_title("h")
cbar = fig.colorbar(cs)
fig.tight_layout()


f, (ax1, ax2) = plt.subplots(1, 2)#generando 2 graficas
im=ax1.imshow(calor)
cs =ax2.contour(contour)
cs1=ax2.contourf(contour)#superponemos dos graficas para tener color y lineas
plt.clabel(cs,inline=1, fontsize=10,colors='black')
for i in range(len(malla)):
    for j in range(len(malla[i])):
        text = ax1.text(j, i, calor[i, j],ha="center",
                       va="center", color="w")
ax1.set_title('Cuadrados')
ax2.set_title('Contorno')
f.tight_layout()
plt.show()