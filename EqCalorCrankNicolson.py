# -*- coding: utf-8 -*-
"""
Created on Wed May  8 12:44:10 2019

@author: Jorge Antonio Matias Lopez
"""

#Ecuación diferencial parcial parabólica(Ecuación de calor), Solución por medio de Crank nicolson
import numpy as np
import math

T = 0.5
l = 1
m = 10
N = 50
k = T/N
h = l/m
#alfa
#lam=(alfa**2)*(k/h**2)
lam = 1

w = np.zeros(m)
lower = np.zeros(m-1)
upper = np.zeros(m-1)
z = np.zeros(m-1)

for i in range(m-1):
    w[i] = math.sin(math.pi*h*(i+1))
print('T=0')
print(w)
lower[0] = 1 + lam
upper[0] = -lam/(2 * lower[0])

for i in range(1,m-2):
    lower[i] = 1 + lam + (lam*upper[i-1])/2
    upper[i]= -lam/( 2 * lower[i])
lower[m-2] = 1 + lam + (lam * upper[m-2])/2

for j in range(1,N+1):
    print(j)
    print(w)
    t = j*k
    z[0] = ((1 - lam) * w[0] + (lam/2) * w[1])/lower[0]
    
    for i in range(1,m-1):
        z[i] = ((1-lam) * w[i] + (lam/2) * (w[i+1] + w[i-1] + z[i-1]))/lower[i]
    w[m-2] = z[m-2]
    
    for i in range(m-3,-1,-1):
        w[i]=z[i]-upper[i]*w[i+1]
    
