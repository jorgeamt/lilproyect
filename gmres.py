# -*- coding: utf-8 -*-
"""
Created on Thu May 30 14:45:18 2019

@author: jor_a
"""
#GMRES
#Datos de entrada=Matriz A(nxn), constante b, 
import numpy as np
import math
def norma(r):#Norma del vector 
    n = r.size
    suma = 0
    for i in range(n):
            #error
            suma = suma + r[i] ** 2
    error=math.sqrt(suma)
    return(error)
#Producto matriz vector
def prodmat(matriz,vector):
    res=[]
    dmat = matriz.shape
    for i in range(dmat[0]):
        suma = 0
        for j in range(dmat[1]):
            suma = suma + matriz[i,j] * vector[j]
        res.append(suma)
    prod=np.array(res)
    return(prod)
###############################################################################
####Producto punto entre dos vectores
def dot(vec1,vec2):
    if vec1.shape == vec2.shape:
        suma = 0
        for i in range(vec1.shape[0]):
            suma = suma + vec1[i] * vec2[i]
        res=np.array(suma)
    return(res)
#Factorización QR
#A partir de a y b obtener c y s
def givens(a,b):
    #print(a,b)
    r=math.sqrt(a ** 2+b ** 2)
    c = a/r
    s = b/r
    return(c,s)
#construir la matriz para realizar la rotación
def matG(i,c,s,m):
    G = np.identity(m)
    j = i
    G[i,j] = c
    G[i-1,j] = -s
    G[i,j-1] = s
    G[i-1,j-1] = c
    return (G)
#Factorización QR 
def QR(R):
    m,n = R.shape
    #print(m,n)
    Q = np.identity(m)
    for j in range(n):
        for i in range(m-1,j,-1):
            if R[i-1,j] != 0 and R[i,j] != 0:
                c,s = givens(R[i-1,j],R[i,j])
                G = matG(i,c,s,m)
                print(G)
                Q = Q.dot(G)
                G = G.transpose()
                #print(R)
                R = G.dot(R)
    return (Q,R)
def gmres(mat,vec):
    x,y = mat.shape
    x0 = np.zeros(x)#vector x inicial
    m = 3
    h = np.zeros([m+1,m])
    V = np.zeros([m,m+1])#Espacio de vectores V
    r = vec-prodmat(mat,x0)#residuo
    normr = norma(r)#Norma del residuo
    print('NORMA',normr)
    v = r/normr#vector unitario v
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
            vj1 = vj1+h[i,j]*vi
        vj1 = -vj1 + Avj
        #Se obtiene V y H 
        h[j+1,j] = norma(vj1)
        V[:,j+1] = vj1/h[j+1,j]

    Q,R = QR(h)#factorización QR de h
    q,q1 = Q.shape
    beta = np.zeros(q)
    beta[0] = normr
    print('beta',beta)
    gk = prodmat(Q,beta)
    print('gk1 \n',gk)
    gk = gk[0:q-1]
    print('gk2 \n',gk)
    Rk = R[0:q-1,:]
    Rinv = np.linalg.inv(Rk)
    print('Rinv \n',Rinv)
    yk = prodmat(Rinv,gk)
    
    for i in range(m):
        x0 = x0 + yk[i]*V[:,i]
    return (V,h,yk,x0,Q,R)
mat = np.array([[2,1,3],
              [3,-1,-2],
              [1,-1,1]],dtype=np.float64)
vector = np.array([7,13,-3],dtype=np.float64)
toler = 0.0000001
V,h,yk,x,Q,R = gmres(mat,vector)
r = vector - prodmat(mat,x)
normr = norma(r)
sol = np.linalg.solve(mat,vector)
#a=np.linalg.inv(A)