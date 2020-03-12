# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 19:06:52 2019

@author: Jorge Antonio Matías López
"""
## Formas
#from math import pi
import numpy as np
import matplotlib.pyplot as plt

# La clase Círculo definida con el método calcArea y el atributo radio
class Circulo:
    
    def __init__(self, radio = 1, centro = (0,0)):
        self.__radio = radio
        self.__centro = centro
        
    
    def setRadio(self, radio):
        self.__radio = radio
        
    def getRadio(self):
        return self.__radio

    def setCentro(self, centro):
        self.__centro = centro
    
    def getCentro(self):
        return self.__centro
        
    def calcArea(self):
        return np.pi * self.__radio ** 2
    
    def dibujar(self):
        p = np.linspace(0, 2 * np.pi, 30)
        x = self.__radio * np.sin(p) + self.__centro[0]
        y = self.__radio * np.cos(p) + self.__centro[1]
        
        fig, ax = plt.subplots(figsize=(8,8))
        ax.plot(x, y,color='blue', lw = 3)
        
        
        
class Rectangulo:
    
    def __init__(self, ladoX = 2, ladoY = 1):
        self.__ladoX = ladoX
        self.__ladoY = ladoY
        
    
    def setladoX(self, ladoX):
        self.__ladoX = ladoX
        
    def getlado(self):
        return self.__ladoX, self.__ladoY

    def setladoY(self, ladoY):
        self.__ladoY = ladoY
    
        
    def calcArea(self):
        
        return self.__ladoX * self.__ladoY
    
    def dibujar(self):
        x = [0,self.__ladoX, self.__ladoX, 0, 0]
        y = [0, 0, self.__ladoY, self.__ladoY, 0]
        
        fig, ax = plt.subplots(figsize=(8,8))
        ax.plot(x, y,color='blue', lw = 3)
        


class Particula:
    
    def __init__(self, posX = 2, posY = 1):
        
        self.__posX = posX
        self.__posY = posY
        
    def setposicion(self, x, y):
        
        self.__posX = x
        self.__posY = y
    
    def getposicion(self):
        
        return self.__posX, self.__posY
    
    
class CampoVelocidad:
       
    def __init__(self, forma = None):
        
        self.forma = forma
        self.lx, self.ly = self.forma.getlado()
        self.nx, self.ny = 20, 20
    
    

    def vel_u(self, x, y, dt=1):
        velx = -np.cos(np.pi * y* 0.1) * np.sin(np.pi * x *0.1)*dt
        return velx

    def vel_v(self, x, y, dt=1):
        vely = np.sin(np.pi * y * 0.1) * np.cos(np.pi * x * 0.1)*dt
        return vely
    
    def graficar(self):
        
        self.xg, self.yg =np.meshgrid(np.linspace(0,self.lx,self.nx),np.linspace(0,self.ly,self.ny))
    
        #Función de velocidad
        u = self.vel_u(self.xg, self.yg)
        v = self.vel_v(self.xg, self.yg)
        #print(u,v)
        
        fig = plt.figure()
        ax = fig.add_subplot(111)
        
        ax.set_aspect('equal')
        ax.quiver(self.xg, self.yg, u, v, units='width')
        #ax.quiverkey(Q, 0.85, 0.9, 1, r'$2 \frac{m}{s}$', labelpos='E',
        #              coordinates='figure')
        
        #ax.streamplot(self.xg, self.yg, u, v,  linewidth=0.5, 
        #              cmap=plt.cm.hot_r, density=1, arrowstyle='->', 
        #              arrowsize=1.0)
        
        
        ax.set_title('Campo vectorial $\\vec{u}$')
        ax.set_xlabel('$x$')
        ax.set_ylabel('$y$')
        
        ax.set_ylim([0, self.ly])
        ax.set_xlim([0, self.lx])
        ax.grid()
        
        return fig, ax
    
    
    


class Trayectoria:
    
    def __init__(self, campo = None, particula = None, dt = 1, pasos = 10):
        
        self.campo = campo
        self.particula = particula
        self.dt = dt
        self.pasos = pasos
    
    
    def calcular(self):
        
        x, y  = self.particula.getposicion()  
        
        self.X = []
        self.Y = []
        
        
        for i in range(self.pasos):
            
            self.X.append(x)
            self.Y.append(y)
            
            x += self.campo.vel_u(x,y,self.dt) 
            y += self.campo.vel_v(x,y,self.dt)
        
        return self.X , self.Y 
    
    def dibujar(self):
        
        
        fig, ax = self.campo.graficar()
        ax.plot(self.X, self.Y, 'b-', marker = 'o')
    
    
        
        
    

circ = Circulo(1, (2,3))

#circ.dibujar()

seg = Rectangulo(10,10)
#seg.dibujar()

part = Particula(3,2)
campvel = CampoVelocidad(seg)
datos = campvel.graficar()

tra = Trayectoria(campvel, part, 0.8, 30)
trayectoria = tra.calcular()
tra.dibujar()
