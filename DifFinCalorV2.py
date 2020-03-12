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
def posicion(n,i,j):#posicion de cada elemento en la matriz de incógnitas 
    pos=(n-1)*(i-1)+j-1
    return pos
###############################################################################
def crs(matriz):
    ###CRS##
    #Al utilizar python se considera el indice inicial = 0
    x,y=matriz.shape#se obtienen las dimensiones de la matriz
    val=[]#lista para almacenar los valores no nulos
    col_ind=[]#lista para almacenar los índices de la columna de cada valor
    ren_elem=[]#Lista para almacenar la informacion de los renglones de la matriz
    cont=0#indice inicial del vector de renglones
    for i in range(x):
        ren_elem.append(cont)
        for j in range(y):
                if matriz[i,j]!=0:
                    val.append(matriz[i,j])#se alnacenan los valores de la matriz
                    col_ind.append(j)#se almacenan los valores de la columna
                    cont=cont+1#incremento de los indices de renglones
    #se transforman las listas a arreglos de numpy
    valor=np.array(val)
    col=np.array(col_ind)
    ren=np.array(ren_elem)
    return(valor,col,ren)
    
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
    if vec1.shape==vec2.shape:#los vectores deben tener el mismo numero de elementos
        suma=0
        for i in range(vec1.shape[0]):
            suma=suma+vec1[i]*vec2[i]
        res=np.array(suma)
    return(res)
#########Modulo de gradiente conjugado con CRS#################################    
def gradc(matriz,vector,tol):
    dim=matriz.shape[0]#tamaño de la matriz
    #Datos a CRS
    val,col_ind,ren_in=crs(matriz)
    #vectores iniciales
    x=np.ones(dim)#vector inicial con valores =1
    r=np.copy(vector)-prodmat(val,col_ind,ren_in,x)#residuo
    d=np.copy(r)#direccion de descenso"d"
    beta=dot(r,r)
    maxitera=1000
    i=0
    t_ini=time()
    while i<maxitera:
        print(i)
        z=prodmat(val,col_ind,ren_in,d)
        rho=beta/(dot(d,z))
        x=x+rho*d
        r=r-rho*z
        gamma=beta
        beta=dot(r,r)
        alfa=(beta/gamma)
        d=r+alfa*d
        #comprobacion
        if np.amax(r)<=tol:
        #if norma(r)<=tol:
            break
        i=i+1
        #print('iteracion',i)
    t_fin=time()
    t_ejecucion=round(t_fin-t_ini,10)
    #print(x)
    print('El tiempo de ejecución Gradiente conjugado es '+str(t_ejecucion)+' segundos')    
    return(x)
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
    #print(res)
    print('El tiempo de ejecución Gauss-Seidel es '+str(t_ejecucion)+' segundos')
    return(res)
############################################################################### 
################################################################################
def gauss_crs(matriz,vector,tol):#Gauss-Seidel usando CRS
    ###CRS
    x,y=matriz.shape
    val=[]#lista para almacenar los valores no nulos
    col_ind=[]#Lista para almacenar el índice de columna de cada dato
    ren_elem=[]#Lista para almacenar los datos de los renglones
    di=[]#lista para almacenar la diagonal
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
                di.append(matriz[i,j])#se extrae la diagonal
    #Transforman las listas en arreglos de numpy           
    valor=np.array(val)
    col=np.array(col_ind)
    ren=np.array(ren_elem)
    diag=np.array(di)
   ##GaussSeidel 
    maxitera=100
    res=np.zeros(x)#arreglo inicial de resultados
    exa=np.linalg.solve(matriz,vector)
    error=[]#lista para evaluar el error
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
    
        del error[:]#Elimina los valores de la lista en cada iteracion
        
            #Comprobación
        for i in range(0,res.size):
            dif=abs(exa[i]-res[i])
            
            error.append(dif)#vector de errores
        if all( i<=tol for i in error) == True:
            break
        #print('Iteracion',k)
    t_fin=time()
    t_ejecucion=round(t_fin-t_ini,10)
    print('El tiempo de ejecución CRS_gauss es '+str(t_ejecucion)+' segundos')
    return(res)#regresa el valor del resultado
def matbld():#construccion de la matriz cuadrada de tamaño n y el vector de constantes
    #en cada ciclo crea un renglón de la matriz y un valor para el vector de constantes
    const=[]#vector de constantes
    mat=[]#renglón de la matriz(se anexa un renglón en cada ciclo)
    for k in range(0,(n-1)**2):#llena el renglón auxiliar con ceros
        mat.append(0)
    pos=posicion(n,i,j)
    mat[pos]=4
    #condiciones de frontera
    
    if j-1>0:
        mat[posicion(n,i,j-1)]=-1#agrega un 1 a la matriz de variables
    else:#Frontera izquierda de f(0,y)=0
        const.append(0)#agrega un 0 al vector de constantes
    if j+1<n:
        mat[posicion(n,i,j+1)]=-1#agrega un 1 a la matriz de variables
    else:#Frontera derecha de f(x,0.5)=200x
        const.append((n-i)*grad)
    
    if i-1>0:
        mat[posicion(n,i-1,j)]=-1
    else:#Frontera superior de f(x,0.5)=200x
        const.append(j*grad)
    
    if i+1<n:
        mat[posicion(n,i+1,j)]=-1
    else:#Frontera inferior de f(x,0)=0
        const.append(0)#agrega un 0 al vector de constantes
    
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
            #err[i,j]=err[i,j]**2
    error=np.amax(err)#norma maximo error
    #error=math.sqrt(np.sum(err))#norma distancia entre vectores
    
    return(error)
###########################################################################################################
###########################################################################################################
####  INICIO  ###  INICIO  ###  INICIO  ###  INICIO  ###  INICIO  ###  INICIO  ###  INICIO  ###  INICIO  ###
e=[]#lista para acumular los errores
tam=[]#Lista de tamaños para n
for k in range(4,9,4):#ciclos de n(tamaño inicial,tamaño final +1,intervalo)
    matriz=[]#lista para generar la matriz de incognitas
    vector=[]#lista para generar el vector de constantes
    mx=100#condición de frontera (valor máximo) 
    n=k#tamaño de la malla considerando el inicio en 0 y con fronteras Dirichlet
    #Las incognitas serán (n-1)^2
    print('n='+str(n))
    tol=0.000000001#tolerancia de solucion
    grad=mx/n#gradiente en las fronteras 
    #Tamaño del dominio
    lx=0.5
    ly=0.5

    dx=lx/n
    dy=ly/n
    #ciclo para construir la matriz de incognitas y el vector de constantes
    for i in range(1,n):
        for j in range(1,n):
            matbld()
    #se transforman las listas a arreglos de numpy  
    mat=np.array(matriz)
    vec=np.array(vector)
    val,col_ind,ren_in=crs(mat)
##############################################################################    
#####soluciones por diferentes metodos########################################
    #solucion=np.linalg.solve(mat,vec)#solución del sistema de ecuaciones linalg python
    #solucion=gradc(mat,vec,tol)#gradiente conjugado con CRS
    #solucion=gauss(mat,vec,tol)#solucion Gauss-Seidel
    solucion=gauss_crs(mat,vec,tol)#solucion Gauss-Seidel matriz CRS
    
    #Arreglos para ordenar las soluciones como una malla
    solcal=np.zeros([(n-1),(n-1)])
    sol=np.zeros([(n-1),(n-1)])
    solex=np.zeros([(n-1),(n-1)])
        
    for i in range(n-1):
        for j in range(n-1):
        
            solex[i,j]=solexa(n,i,j)#matriz solución exacta
            solcal[i,j]=solucion[(i*(n-1))+j]#matriz solucion calculada
            sol[i,j]=round(solucion[(i*(n-1))+j],2)#matriz solución calculada con redondeo para graficar
    errorN2=error(solex,solcal,n)#cálculo del error
###############################################################################
#### Cálculo de la convergencia###############################################
    #e.append(errorN2)
    #tam.append(n)
    
    #e.append(math.log10(errorN2))#Log-Log
    e.append(-math.log10(errorN2))#-Log-Log
    tam.append(math.log10(n))

m=[]#pendientes
for i in range(len(e)-1):
    d_e=e[i]-e[i+1]
    d_tam=tam[i]-tam[i+1]
    m.append(d_e/d_tam)

X1=np.array(tam)
Y1=np.array(e)

#razón de convergencia
f, (ax1) = plt.subplots()
ax1.plot(X1,Y1)
ax1.grid()
ax1.set(xlabel='tamaño', ylabel='Error',
       title='Razón de convergencia')

###############################################################################
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