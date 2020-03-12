# -*- coding: utf-8 -*-
"""
Created on Tue May 14 13:31:41 2019

@author: jor_a
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import math

def crs(diagP,diagS,dim):
    v=[diagP,diagS]
    inter=[diagS,diagP,diagS]
    finv=[diagS,diagP]
    col=[0,1]
    finc=[dim-2,dim-1]
    cont=0
    for i in range(n-2):
        v=v+inter
        for j in range(3):
            col.append(cont+j)
        cont=cont+1
    v=v+finv
    col=col+finc
    ren_in=np.zeros(dim,dtype=int)
    ren_in[0]=0
    ren_in[1]=2
    for i in range(dim-2):
        ren_in[i+2]=ren_in[i+1]+3
    valor=np.array(v)
    col_ind=np.array(col)
    return(valor,col_ind,ren_in) 
    

def prodmat(valor,col,ren,vector):##Producto matricial utilizando CRS
    res=[]#lista para almacenar los resultados del producto
    
    for i in range(0,ren.size):
        suma=0
        if i != ren.size-1:
            for j in range(ren[i],ren[i+1]):
                suma=suma+valor[j]*vector[col[j]]
            res.append(suma)
        else:
            for j in range(ren[i],valor.size):
                suma=suma+valor[j]*vector[col[j]]               
            res.append(suma)
    prod=np.array(res)
    return(prod) 
    
def dot(vec1,vec2):
    if vec1.shape==vec2.shape:#los vectores deben tener el mismo numero de elementos
        suma=0
        for i in range(vec1.shape[0]):
            suma=suma+vec1[i]*vec2[i]
        res=np.array(suma)
    return(res)
    
def gradc(val,col_ind,ren_in,vector,tol):#datos de entrada (matriz CRS, vector, tolerancia)
    dim=vector.size#tamaño del vector solucion
    #vectores iniciales
    x=np.ones(dim)#vector inicial con valores =1
    r=np.copy(vector)-prodmat(val,col_ind,ren_in,x)#residuo
    d=np.copy(r)#direccion de descenso"d"
    beta=dot(r,r)
    maxitera=100000
    i=0
    while i<maxitera:
        z=prodmat(val,col_ind,ren_in,d)
        rho=beta/(dot(d,z))
        x=x+rho*d
        r=r-rho*z
        gamma=beta
        beta=dot(r,r)
        alfa=(beta/gamma)
        d=r+alfa*d
        #comprobacion
        if np.amax(r) <= tol:#norma del error maximo
            break
        i=i+1
        #print('iteracion',i)    
    return(x)
###############################################################################
    
T = 0.5 #Tiempo máximo
l = 1 #longitud del dominio
m = 10 #Número de nodos
N = 50 #pasos de tiempo
k = T/N #dt
h = l/m #dx
#alfa #
#lam=(alfa**2)*(k/h**2)
lam = 1
tol = 0.0000001
n = m-1#número de incógnitas

diagA = 1+lam
diagB = 1-lam
a = -lam/2
b = -a
valA,col_indA,ren_inA = crs(diagA,a,n)
valB,col_indB,ren_inB = crs(diagB,b,n)


w = np.zeros(n)
res = np.zeros(n+2)
x = np.arange(0, l+h, h)
y=[]
for i in range(n):
    w[i] = math.sin(math.pi*h*(i+1))
res[1:-1] = w[:]
#print(w)

for j in range(1,N+1):
    
    t = k*j
    vec = prodmat(valB,col_indB,ren_inB,w)
    w = gradc(valA,col_indA,ren_inA,vec,tol)
    res[1:-1] = w[:]
    #print(res)
    y.append(list(res[:]))
    #print(y[-1])
    #print('t=',t)
    
    
    
#imprimir grafica con animación
fig, ax = plt.subplots()
xdata, ydata = [x[:]] , [res[:]]
ax.set(xlabel = 'Dominio x', ylabel = 'Temperatura',
       title='Difusión de calor') 
ln, = plt.plot([], [])


def init():
    ax.set_xlim(0, l+h)
    ax.set_ylim(0, 1)
    return ln,

def update(paso):
    print(paso)
    xdata = x[:]  
    ydata = paso[:]
    ln.set_data(xdata, ydata)
    return ln,

ani = FuncAnimation(fig, update, y ,
                    init_func=init, blit=True, repeat= False)
plt.show()





