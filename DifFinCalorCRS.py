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
def matbld(n):#construccion de la matriz cuadrada de tamaño n a CRS, y el vector de constantes 
    valor=[]#lista de valores no nulos de la matriz
    col_ind=[]#Lista de indices de las columnas de los valores
    ren=[]#Lista de Elementos en cada renglón
    vector=[]#Construcción del vector de valores
    cont=0#Elemento inicial de la lista de elementos de renglón
    #se recorre la matriz dentro de las condiciones de frontera
    for i in range(1,n):#0<i<n
        for j in range(1,n):#0<j<n
            ren.append(cont)#se agrega el elemento a la lista de elementos del renglón
            const=[]#vector de constantes
            col=[]#lista auxiliar para agregar los índices de columnas
            val=[]#Lista auxiliar para agregar los valores no nulos 
            pos=posicion(n,i,j)#indice de posición de la incógnita en i,j
            col.append(pos)#Se agrega el indice de la posicion a la lista
            val.append(4)#En la posición i,j el valor es 4, y se agrega
            #condiciones de frontera
            
            if j-1>0:#condicion de frontera izquierda 
                #si no está en contacto con la frontera
                if posicion(n,i,j-1)>pos:#si el indice es mayor que el de i,j
                    
                    col.append(posicion(n,i,j-1))#se agrega la posición al final
                    val.append(-1)#Se agrega al final -1 en la lista auxiliar
                else:#si el indice es menor que el de i,j
                    col.insert(0,posicion(n,i,j-1))#se agreja la posición al inicio
                    val.insert(0,-1)#se agrega el valor -1 al inicio en la lista auxiliar 
            else:#si está en contacto con la frontera izquierda
                const.append(0)#Se agrega un 0 a la lista de constantes
            if j+1<n:#condicion de frontera derecha
                #si no está en contacto con la frontera
                if posicion(n,i,j+1)>pos:#si el indice es mayor que el de i,j
                    col.append(posicion(n,i,j+1))#se agrega la posición al final
                    val.append(-1)#Se agrega al final -1 en la lista auxiliar
                else:#si el indice es menor que el de i,j
                    col.insert(0,posicion(n,i,j+1))#se agreja la posición al inicio
                    val.insert(0,-1)#se agrega el valor -1 al inicio en la lista auxiliar 
            else:#si está en contacto con la frontera derecha
                const.append((n-i)*grad)##se agrega el valor de la frontera en ese renglón en la lista constantes
            if i-1>0:#condicion de frontera superior
                #si no está en contacto con la frontera
                if posicion(n,i-1,j)>pos:#si el indice es mayor que el de i,j
                    col.append(posicion(n,i-1,j))#se agrega la posición al final
                    val.append(-1)#Se agrega al final -1 en la lista auxiliar
                else:#si el indice es menor que el de i,j
                    col.insert(0,posicion(n,i-1,j))#se agreja la posición al inicio
                    val.insert(0,-1)#se agrega el valor -1 al inicio en la lista auxiliar 
            else:#si está en contacto con la frontera superior
                const.append(j*grad)#se agrega el valor de la frontera en esa columna en la lista constantes
            if i+1<n:#condicion de frontera inferior
                #si no está en contacto con la frontera
                if posicion(n,i+1,j)>pos:#si el indice es mayor que el de i,j
                    col.append(posicion(n,i+1,j))#se agrega la posición al final
                    val.append(-1)#Se agrega al final -1 en la lista auxiliar
                else:#si el indice es menor que el de i,j
                    col.insert(0,posicion(n,i+1,j))#se agreja la posición al inicio
                    val.insert(0,-1)#se agrega el valor -1 al inicio en la lista auxiliar 
            else:#si está en contacto con la frontera inferior
                const.append(0)#se agrega el valor 0 en la lista constantes
            col_ind=col_ind+col#se agregan los indices de columnas a la lista principal
            valor=valor+val#se agregan los valores a la lista de valores principal
            cont=cont+len(val)#se suma a cont la cantidad de elementos de la lista auxiliar val
            vector.append(sumalista(const))#se agrega la suma de las constantes obtenidas al vector de constantes
    #Se transforman las listas a arreglos de numpy 
    valr=np.array(valor)
    col_inx=np.array(col_ind)
    ren_in=np.array(ren)
    vec=np.array(vector)
    return(valr,col_inx,ren_in,vec)

###############################################################################    
def prodmat(valor,col,ren,vector):##Producto matricial utilizando CRS
    res=[]#lista para almacenar los resultados del producto
    
    for i in range(0,ren.size):#0<i<ren
        suma=0#suma de los productos 
        if i != ren.size-1:#para evitar excepciones en la evaluación se usan los elementos de ren 0-ren-2
            for j in range(ren[i],ren[i+1]):
                suma=suma+valor[j]*vector[col[j]]#se suman los productos vector-columna
            res.append(suma)#Se agrega el valor final de la suma a la lista de resultados
        else:#si se usa el indice ren-1 se utiliza la información del número de elementos de val
            for j in range(ren[i],valor.size):
                suma=suma+valor[j]*vector[col[j]]#se suman los productos vector-columna       
            res.append(suma)#Se agrega el valor final de la suma a la lista de resultados
    prod=np.array(res)#se transforma en un arreglo de numpy
    return(prod) 
####Producto punto entre dos vectores
def dot(vec1,vec2):
    if vec1.shape==vec2.shape:#los vectores deben tener el mismo numero de elementos
        suma=0
        for i in range(vec1.shape[0]):
            suma=suma+vec1[i]*vec2[i]
        res=np.array(suma)
    return(res)
##############################################################################
def solexa(n):#solucion exacta
    #Tamaño del dominio
    lx=0.5
    ly=0.5

    dx=lx/n
    dy=ly/n
    solex=[]#lista para almacenar los resultados de la solución exacta
    for i in range(n-1):
        for j in range(n-1):
            x=(j+1)*dx#se calcula x & su valor en el plano cartesiano
            y=(n-i-1)*dy#se calcula y & su valor en el plano cartesiano
            exa=400*x*y
            solex.append(exa)
    exacta=np.array(solex)#Se convierte en un arreglo de numpy
    
    return(exacta)
##########################################################################################################
def error(exa,cal):
    #cálculo del error
    err=np.zeros(exa.size)#vector de errores lleno de ceros
    for i in range(n-1):
        err[i]=abs(exa[i]-cal[i])#error
            #err[i,j]=err[i,j]**2
    error=np.amax(err)#norma maximo error
    #error=math.sqrt(np.sum(err))#norma distancia entre vectores
    return(error)
################################################################################
def malla(vector):
    #Arreglos para ordenar las soluciones como una malla
    n=int(math.sqrt(vector.size))
    malla=np.zeros([n,n])
    for i in range(n):
        for j in range(n):
            malla[i,j]=round(vector[(i*(n))+j],2)#matriz solucion calculada
            #malla[i,j]=vector[(i*(n))+j]#matriz solucion calculada
    return(malla)
#########Modulo de gradiente conjugado con CRS#################################    
def gradc(val,col_ind,ren_in,vector,tol):#datos de entrada (matriz CRS, vector, tolerancia)
    dim=vector.size#tamaño del vector solucion
    n=int(math.sqrt(dim)+1)
    exa=solexa(n)
    #vectores iniciales
    x=np.ones(dim)#vector inicial con valores =1
    r=np.copy(vector)-prodmat(val,col_ind,ren_in,x)#residuo
    d=np.copy(r)#direccion de descenso"d"
    beta=dot(r,r)
    maxitera=100000
    i=0
    t_ini=time()
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
        err=error(exa,x)
        if err<=tol:
        #if np.amax(r)<=tol:#norma del error maximo
            break
        i=i+1
        #print('iteracion',i)
    t_fin=time()
    t_ejecucion=round(t_fin-t_ini,10)
    #print(x)
    print('El tiempo de ejecución Gradiente conjugado es '+str(t_ejecucion)+' segundos')    
    return(x)
##################################################################### 
################################################################################

   ##GaussSeidel 
def gauss_crs(valor,col,ren,vector,tol):
    n=int(math.sqrt(vector.size)+1)
    exa=solexa(n)
    ##GaussSeidel 
    maxitera=100000
    res=np.zeros(vector.size)#arreglo inicial de resultados
    k=0
    t_ini=time()
    while k<maxitera:
        suma=0
        k=k+1
        d=0
        for i in range(0,ren.size):
            suma=0
            if i != ren.size-1:
                for j in range(ren[i],ren[i+1]):
                    if col[j]!=d:
                        suma=suma+valor[j]*res[col[j]]   
                    else:
                        diag=valor[j]
                res[i]=(vector[i]-suma)/diag
            else:
                for j in range(ren[i],valor.size):
                    if col[j]!=d:
                        suma=suma+valor[j]*res[col[j]] 
                    else:
                        diag=valor[j]
                res[i]=(vector[i]-suma)/diag
            d=d+1
        err=error(exa,res)
        if err<=tol:
            break
    t_fin=time()
    t_ejecucion=round(t_fin-t_ini,10)
    print('El tiempo de ejecución CRS_gauss es '+str(t_ejecucion)+' segundos')
    return(res)

###########################################################################################################
###########################################################################################################
####  INICIO  ###  INICIO  ###  INICIO  ###  INICIO  ###  INICIO  ###  INICIO  ###  INICIO  ###  INICIO  ###
e=[]#lista para acumular los errores
e1=[]
tam=[]#Lista de tamaños para n
for k in range(31,32,3):#ciclos de n(tamaño inicial,tamaño final +1,intervalo)
    mx=100#condición de frontera (valor máximo) 
    n=k#tamaño de la malla considerando el inicio en 0 y con fronteras Dirichlet
    #Las incognitas serán (n-1)^2
    print('n='+str(n),'incógnitas='+str((n-1)**2) )
    tol=0.000000001#tolerancia de solucion
    grad=mx/n#gradiente en las fronteras 
    valor,col_in,ren_in,vec=matbld(n)
    
##############################################################################    
#####soluciones por diferentes metodos########################################
    solucion=gradc(valor,col_in,ren_in,vec,tol)#gradiente conjugado con CRS
    #solucion=gauss_crs(valor,col_in,ren_in,vec,tol)#solucion Gauss-Seidel matriz CRS
    
    solex=solexa(n)
    errorN2=error(solex,solucion)#cálculo del error
            #sol[i,j]=round(solucion[(i*(n-1))+j],2)#matriz solución calculada con redondeo para graficar
    
###############################################################################
#### Cálculo de la convergencia###############################################
#    e1.append(errorN2)
#    tam.append(n)
    
    #e.append(math.log10(errorN2))#Log-Log
    e.append(-math.log10(errorN2))#-Log-Log
    tam.append(math.log10(n))

m=[]#pendientes
for i in range(len(e)-1):
    d_e = e[i]-e[i+1]
    d_tam = tam[i]-tam[i+1]
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
sol = malla(solucion)
#graficas
fig, ax = plt.subplots()#Solo genera una grafica
im = ax.imshow(sol)
'''for i in range(sol.shape[0]):
    for j in range(sol.shape[0]):
        text = ax.text(j, i, sol[i, j],ha="center",
                       va="center", color="w")'''
ax.set_title("Tamaño de malla=31")
cbar = fig.colorbar(im)
fig.tight_layout()
plt.show