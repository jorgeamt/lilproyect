# -*- coding: utf-8 -*-
"""
Created on Sat Feb  2 10:24:15 2019

@author: Jorge Antonio Matias López
"""
#Programa que soluciona la ecuación de laplace mediante diferencias finitas, ejemplo del libro Numerical Analysis
#capitulo 12.1 generalizado para n dimensiones
import numpy as np
import matplotlib.pyplot as plt
import math
from time import time

####################################################################

def sumalista(listaNumeros):#suma los elementos de una lista
    Suma = 0
    for i in listaNumeros:
        Suma = Suma + i
    return Suma
####################################################################
def posicion(n,i,j):#posicion en la matriz
    pos=(n-1)*(i-1)+j-1
    return pos
#####################################################################
def gauss(matriz,vect,tol):#Método de Gauss-Seidel
    print('Gauss-Seidel')
    dim=matriz.shape
    comp=np.zeros(dim[0])
    itera=1000
    res=np.zeros(dim[0])
    error=[]
    k=0
    t_ini=time()
    while k<itera:
        suma=0
        k=k+1
        for ren in range(0,dim[0]):
            suma=0
            for col in range(0,dim[1]):
                if (col != ren):
                    suma=suma+matriz[ren,col]*res[col]               
            res[ren]=(vect[ren]-suma)/matriz[ren,ren]
        del error[:]
            #Comprobación
        for ren in range(0,dim[0]):
            suma=0
            for col in range(0,dim[1]):
                suma=suma+matriz[ren,col]*res[col]  
            comp[ren]=suma
            dif=abs(comp[ren]-vect[ren])
            
            error.append(dif)

        #print('Iteracion',k)
        if all( i<=tol for i in error) == True:
            break
    t_fin=time()
    t_ejecucion=round(t_fin-t_ini,10)
    print(res)
    print('El tiempo de ejecución Gauss-Seidel es '+str(t_ejecucion)+' segundos')
    return(res)
############################################################################### 
def gauss_crs(matriz,vector,tol):#Gauss-seidel con CRS
    print('CRS Gauss-Seidel')
    x,y=matriz.shape
    val=[]
    col_ind=[]
    ren_elem=[]
    di=[]
    cont=0
    for i in range(x):
        ren_elem.append(cont)
        for j in range(y):
            if i!=j:
                if matriz[i,j]!=0:
                    val.append(matriz[i,j])
                    col_ind.append(j)
                    cont=cont+1
            else:
                di.append(matriz[i,j])
                
    valor=np.array(val)
    col=np.array(col_ind)
    ren=np.array(ren_elem)
    diag=np.array(di)
    #return(valor,col,ren,diag) 
   ##GaussSeidel 
    maxitera=1000
    res=np.zeros(x)
    exa=np.linalg.solve(matriz,vector)
    error=[]
    k=0
    t_ini=time()
    while k<maxitera:
        suma=0
        k=k+1
        for i in range(0,ren.size):
            suma=0
            if i != ren.size-1:
                for j in range(ren[i],ren[i+1]):
                    suma=suma+valor[j]*res[col[j]]  
                    
                res[i]=(vector[i]-suma)/diag[i]
            else:
                for j in range(ren[i],valor.size):
                    suma=suma+valor[j]*res[col[j]]               
                res[i]=(vector[i]-suma)/diag[i]
    
        del error[:]
        
        
            #Comprobación
        for i in range(0,res.size):
            dif=abs(exa[i]-res[i])
            error.append(dif)
        if all( i<=tol for i in error) == True:
            break
    t_fin=time()
    t_ejecucion=round(t_fin-t_ini,10)
    #print(res)
    print('El tiempo de ejecución CRS_gauss es '+str(t_ejecucion)+' segundos')
    return res
################################################################################
def matbld():#construccion de la matriz cuadrada de tamaño n y el vector de constantes
    const=[]#vector de constantes
    mat=[]#renglón de la matriz(se anexa un renglón en cada ciclo)
    for k in range(0,(n-1)**2):#llena el renglón auxiliar con ceros
        mat.append(0)
    pos=posicion(n,i,j)
    mat[pos]=4
    #condiciones de frontera
    if j-1>0:
        mat[posicion(n,i,j-1)]=-1
    else:
        const.append(0)
    if j+1<n:
        mat[posicion(n,i,j+1)]=-1
    else:
        const.append((n-i)*grad)
    if i-1>0:
        mat[posicion(n,i-1,j)]=-1
    else:
        const.append(j*grad)
    if i+1<n:
        mat[posicion(n,i+1,j)]=-1
    else:
        const.append(0)
    
    matriz.append(mat)#agregla el renglón creado a la matriz principal
    #agrega la suma de los coeficientes en las fronteras al vector de coeficientes
    vector.append(sumalista(const))
    return()
##########################################################################################################
def solexa(n,i,j):#solucion exacta
    x=(j+1)*dx
    y=(n-i-1)*dy
    exa=abs(400*x*y)
    '''print('x','y','Solución exacta')
    print(x,y,exa)'''
    return(exa)
##########################################################################################################
def error(exa,cal,n):
    err=np.zeros([(n-1),(n-1)])
    for i in range(n-1):
        for j in range(n-1):
            err[i,j]=abs(exa[i,j]-cal[i,j])#error
            err[i,j]=err[i,j]**2
    error=math.sqrt(np.sum(err))
    
    return(error)
###########################################################################################################
###########################################################################################################
####  Inicio  ###
e=[]
tam=[]
for k in range(4,17,2):
    matriz=[]
    vector=[]
    mx=100#condición de frontera (valor máximo) 
    n=k
    print('n='+str(n))
    tol=0.00000001
    grad=mx/n
    #Tamaño del dominio
    lx=0.5
    ly=0.5

    dx=lx/n
    dy=ly/n
    for i in range(1,n):
        for j in range(1,n):
            matbld()
        
    mat=np.array(matriz)
    vec=np.array(vector)
    #soluciones por diferentes metodos
    solucion=np.linalg.solve(mat,vec)#solución del sistema de ecuaciones linalg python
    #sol_gauss=gauss(mat,vec,tol)#solucion Gauss-Seidel
#    solucion=gauss_crs(mat,vec,tol)#solucion Gauss-Seidel matriz CRS

    solcal=np.zeros([(n-1),(n-1)])
    sol=np.zeros([(n-1),(n-1)])

    solex=np.zeros([(n-1),(n-1)])


    
    for i in range(n-1):
        for j in range(n-1):
        
            solex[i,j]=solexa(n,i,j)#matriz solución exacta
            solcal[i,j]=solucion[(i*(n-1))+j]#matriz solucion calculada
            sol[i,j]=round(solucion[(i*(n-1))+j],2)#matriz solución calculada
        

    
    errorN2=error(solex,solcal,n)

#    e.append(errorN2)
#    tam.append(n)
    
    e.append(math.log10(errorN2))
    tam.append(math.log10(n))
    
#    e.append(-math.log10(errorN2))
#    tam.append(math.log10(n))
    X=np.array(tam)
    Y=np.array(e)
#razón de convergencia
fig,ax0=plt.subplots()
ax0.plot(X,Y)
ax0.grid()
ax0.set(xlabel='tamaño', ylabel='Error',
       title='Razón de convergencia -log-log')

#graficas
'''fig, ax = plt.subplots()#Solo genera una grafica
im = ax.imshow(sol)
for i in range(sol.shape[0]):
    for j in range(sol.shape[0]):
        text = ax.text(j, i, sol[i, j],ha="center",
                       va="center", color="w")
ax.set_title("Solución °C")
cbar = fig.colorbar(im)
fig.tight_layout()'''
plt.show