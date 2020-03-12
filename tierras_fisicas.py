# -*- coding: utf-8 -*-
"""
Created on Fri Sep 27 09:22:18 2019

@author: Jorge Antonio Matías López
"""

'''
Programa que a partir de datos de SEV realiza los cálculos para determinar la Resistencia mínima 
para hacer un arreglo de tierras físicas.
'''
import numpy as np
import math
print('''
Cálculo de tierras físicas
      ''')
dat = input('¿Desea importar los datos desde un archivo de texto o manualmente? m=manual, cualquier otra tecla= importar \n')

if dat == 'm':
    capas = int(input('ingrese el número de capas: \n'))
    E = []
    rho = []
    
    for i in range(capas):
        print(f'Valores de la Capa {i+1}')
        espesor = float(input('Ingrese el valor del espesor en m \n'))
        E.append(espesor)
        resistividad = float(input('Ingrese el valor de resistividad en ohm-m \n'))
        rho.append(resistividad)
    E = np.array(E)
    rho = np.array(rho)
else:
    file = input(
'''
Datos desde archivo: 
Introduzca el nombre del archivo con extensión(txt):
Formato:
Espesor    Resistividad (Encabezado)
E0         rho0
E1         rho1
 .            .
 .            .
 \n
''')
    datos = np.loadtxt(file, skiprows=1, dtype=np.float32)
    E = datos[:, 0]
    rho = datos[:, 1]
    capas = datos.shape[0]
T = 0 #Resistencia Transversal unitaria
S = 0 #Condictancia longitudinal

for i in range(capas):
    T += E[i] * rho[i]
    S += E[i] / rho[i]
    
rhoeq = math.sqrt(T/S) #Cálculo de resistividad equivalente

elec = input('¿Desea modificar las características del electrodo?(por defecto long=3.0m y radio=0.0079m) Y/N \n')
if elec == 'Y':
    a = float(input('Introduce el radio del electrodo en m: \n'))
    l = float(input('Introduce la longitud del electrodo en m: \n'))
else:
    a = 0.0079
    l = 3.0
    
R = rhoeq/(2*math.pi*l)* math.log(2*l/a)#resistencia

if R <= 10:
    print(f'La resistencia calculada para un electrodo es {R} y está debajo del umbral de 10 ohm')
    print('¡Adios!')
else:
    print(f'La resistencia obtenida ha sido {R}')
    elec = input('¿Desea calcular un arreglo en paralelo? Y/N \n')
    if elec == 'Y':
        d = float(input('Ingrese la distancia entre electrodos en m: \n'))
        n = 2
        while True:
            A = math.pow(a*math.pow(d,n-1),1/n)#Cálculo del radio equivalente para cada iteración
            R = rhoeq/(2*math.pi*l)* math.log(2*l/A)
            
            if R <= 10:
                print(f'La resistencia calculada para un arreglo en paralelo con {n} electrodos es {R}\ny está debajo del umbral de 10 ohm')
                print('¡Adios!')
                break
            elif n >= 5:
                print(f'Se ha alcanzado el número de electrodos disponibles, la resistencia obtenida es {R}')
                print('¡Adios!')
                break
            n += 1
        
    else:
        print('Fin del programa, ¡Adios!')
