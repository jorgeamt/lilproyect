# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 23:55:34 2019

@author: jor_a
"""
dx=1
dy=1
dt=1
S=0.4
I=0.001
h0=1
T=0.1
malla=[100, 100, 100, 100],[h0,h0,h0,h0],[h0,h0,h0,h0],[0,h0,h0,h0]
alfa=(T*dt)/(S*dx*dx)
beta=(T*dt)/(S*dy*dy)
simp=(I*dt)/S
#h=malla[i][j]+simp+alfa*(malla[i][j-1]+malla[i][j+1]-2*malla[i][j])+beta*(malla[i-1][j]+malla[i+1][j]-2*malla[i][j])

#iteracion=input('Número de iteraciones')
ciclos=60
cont=1
while cont<=ciclos:
      
    for i in range(len(malla)):
        for j in range(len(malla[i])):
            
            if len(malla)-1>i>0:#condición de frontera S1
                if j==0:#C. F. izquierda
                    malla[i][j]=round((malla[i][j]+simp+alfa*(2*malla[i][j+1]-2*malla[i][j])+beta*(malla[i-1][j]+malla[i+1][j]-2*malla[i][j])),5)
                    
                elif j==len(malla[i])-1:#C.F. derecha
                    malla[i][j]=round((malla[i][j]+simp+alfa*(2*malla[i][j-1]-2*malla[i][j])+beta*(malla[i-1][j]+malla[i+1][j]-2*malla[i][j])),5)
                    
                else:
                    malla[i][j]=round((malla[i][j]+simp+alfa*(malla[i][j-1]+malla[i][j+1]-2*malla[i][j])+beta*(malla[i-1][j]+malla[i+1][j]-2*malla[i][j])),5)
                    
            elif i==len(malla)-1:#frontera inferior
                if j==0:#C. F. izquierda
                    pass#condicion de frontera S2
                      
                elif j==len(malla[i])-1:#C.F. derecha
                    malla[i][j]=round((malla[i][j]+simp+alfa*(2*malla[i][j-1]-2*malla[i][j])+beta*(2*malla[i-1][j]-2*malla[i][j])),5)
                    cont=cont+1
                    
                else:
                    malla[i][j]=round((malla[i][j]+simp+alfa*(malla[i][j-1]+malla[i][j+1]-2*malla[i][j])+beta*(2*malla[i-1][j]-2*malla[i][j])),5)

