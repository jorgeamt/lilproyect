# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 16:08:56 2019

@author: Jorge Antonio Matías López
"""

'''
Proyecto python STUMOD 

Parámetros hidraulicos

HLR      cm/day      Hydraulic loading rate
α1        adim       Parameter α in Gradner’s analytical equation for pressure distribution (also referred to as αG)
α2        adim       Parameter α in the soil water retention function (also referred to as αVG)
Ks       cm/day      Saturated hydraulic conductivity (also referred to as Ksat)
θ1        adim       Residual soil moisture (also referred to as θr)
θ2        adim       Saturated soil moisture (also referred to as θs)
n         adim       Parameter n in the soil water retention function
M         adim       Parameter m in the soil water retention function
L         adim       Tortuosity parameter

Biomat parameters

Kb       cm/day      Biomat hydraulic conductivity
BT         cm        Biomat thickness

Effluent quality parameters

Co-NH4    mg-N/L     Effluent ammonium-N concentration
Co-NO3    mg-N/L     Effluent nitrate-N concentration

Nitrification parameters

Kr-max  mg-N/L-day   Maximum nitrification rate
Km-nit   mg-N/L      Half-saturation constant for ammonium-N
e2       adim         Empirical exponent for nitrification
e3       adim         Empirical exponent for nitrification
fs       adim         Value of the soil water response function at saturation
fwp      adim         Value of the soil water response function at wilting point
swp      adim         Relative saturation at wilting point
sl       adim         Relative saturation for biological process (lower limit)
sh       adim         Relative saturation for biological process (upper limit)
β1       adim         Empirical coefficient for temperature function for nitrification (also referred to as βnit)
Topt     °C           Optimum temperature for nitrification


Denitrification Parameters


Vmax  mg-N/L-day Maximum denitrification rate
Km-dnt mg-N/L Half-saturation constant for nitrate-N
ednt adim Empirical exponent for denitrification
sdn adim A threshold relative saturation (dimensionless)
β2 adim An empirical coefficient for temperature function(also referred to as βdnt)
Topt °C Optimum temperature for denitrification
αc adim An empirical exponent for carbon content adjustment


Ammonium Sorption Parameters
kd L/kg Adsorption isotherm
Ρ kg/L Soil bulk density

Temperature parameters
T °C Soil temperature

Treatment Depth
D cm Soil depth

Output values at treatment depth
C/Co NH4 mg-N L−1 Fraction of ammonium-N remaining at soil depth D
C/Co TotN mg-N L−1 Fraction of total N remaining at soil depth D

'''
import numpy as np
import matplotlib.pyplot as plt
import math

##Orientado a objetos, el objeto es modelo 
class Modelo:
    
    def __init__(self, D = 100, dz = 0.5, propiedades = None):
        self.D = D
        self.dz = dz    
        self.prop = propiedades
        self.prof = None
        self.y =  None
        self.hum = None
        self.Se = None
        self.K = None
        self.vz = propiedades['HLR']/propiedades['qr']
        self.fz = None
        self.Td = self.dz/self.vz
        self.WFP = None
        self.R = None
        self.krmax = None
        self.Vmax = None
        self.C_NH4 = None
        self.C_NO3 = None
        self.Ctotal = None
        
    
    def z(self):
        D = self.D
        dz = self.dz
        intervalo = int(D//dz)
        prof = np.linspace(0, D, intervalo)
        self.prof = prof
        return prof
    
    
    def psi(self, z):
        prop = self.prop
        aG = prop['aG']
        Yo = prop['Yo']
        Ks = prop['Ks']
        HLR = prop['HLR']
        #y = (1/aG)*math.log(math.exp(aG*(Yo-z)) + (HLR/Ks)*(math.exp(-aG*z)-1))
        log = (HLR/Ks)*math.exp(-aG*Yo)*(math.exp(aG*z)-1)+1
        y = Yo - z + (1/aG) * math.log(log)
        return y
        
    
    def tetha(self):####### intentar con map, tetha(self, y)
        psi = self.y
        aVG = self.prop['aVG']
        qs = self.prop['qs']
        qr = self.prop['qr']
        n = self.prop['n']
        m = 1 - 1/n
        print(aVG, qs, qr,n,m, sep='\n')
        tetha = np.zeros(psi.shape[0])
        
        for ind, val in enumerate(psi):
            if val > 0:
                tetha[ind] = qs
                
                
            else:
                ss = qr + (qs-qr)/ math.pow((1 + math.pow(abs(aVG*val),n)),m)
                tetha[ind] = ss
        self.hum = tetha
        return tetha

    
    def se(self, tetha):
        qr = self.prop['qr']
        qs = self.prop['qs']
        Se = (tetha - qr) / (qs - qr)
        return Se
    
    
    def k(self, se):
        Ks = self.prop['Ks']
        l = self.prop['l']
        n = self.prop['n']
        m = 1 - 1/n
        try:
            
            k = Ks*math.pow(se,l)*math.pow(1-math.pow(1-math.pow(se, 1/m),m),2)
        except:
            k=0
            print(se, m)
        return k
    
    def vz(self):
        pass
    
    
    def fWFP(self):### intentar con map fWFP(self, hum)
        
        qs = self.prop['qs']
        swp = self.prop['swp']
        hum = self.hum
        WFP = np.zeros(len(hum))
        for i, tetha in enumerate(hum):
            sw = tetha/qs
            if sw > swp:
                WFP[i] = sw
            else:
                WFP[i] = swp
        self.WFP = WFP
        return WFP
    
    def fcont_c(self):
        ac =  self.prop['ac']
        cont_c = np.zeros(len(self.prof))
        for i, z in enumerate(self.prof):
            cont_c[i] = math.exp(-ac*z)
        self.fz = cont_c
        return cont_c
    
    
    def fswnt(self):
        WFP = self.WFP
        fs = self.prop['fs']
        fwp = self.prop['fwp']
        swp = self.prop['swp']
        e2 = self.prop['e2']	 # e2
        e3 = self.prop['e3']
        sl = self.prop['sl']
        sh = self.prop['sh']
        
        fsw_nt = np.zeros(len(WFP))
        
        for i, wfp in enumerate(WFP):
            if sh < wfp <= 1:
                fsw_nt[i] = fs + (1-fs)*math.pow((1-wfp)/(1-sh),e2)
            elif sl <= wfp <= sh:
                fsw_nt[i] = 1
            elif swp <= wfp <= sl:
                fsw_nt[i] = fwp + (1-fwp)*math.pow((wfp-swp)/(sl-swp),e3)
        self.fsw_nt = fsw_nt
        return fsw_nt
        
      
    def fswdnt(self):
        WFP = self.WFP
        sdn = self.prop['sdn']
        ednt = self.prop['ednt']
        
        fsw_dnt = np.zeros(len(WFP))
        
        for i, wfp in enumerate(WFP):
            if wfp < sdn:
                fsw_dnt[i] = 0
                
            else:
                fsw_dnt[i] = math.pow((wfp-sdn)/(1-sdn),ednt)
        self.fsw_dnt = fsw_dnt
        return fsw_dnt
    
    
    def fR(self):
        kd = self.prop['kd']
        tetha = self.hum
        rho =  self.prop['rho']
        
        r = np.zeros(len(tetha))
        
        for i, h in enumerate(tetha):
            r[i] = 1 + (rho*kd/h)
        self.R = r
        return r
    
    def f_krmax(self):
        kr_max = self.prop['kr_max']
        kr = self.fsw_nt * kr_max
        self.krmax = kr
        return kr
    
    
    def fVmax(self):
        Vmax = self.prop['Vmax']
        Vm = self.fsw_dnt * Vmax
        self.Vmax = Vm
        return Vm
    
    def fconc(self, km, Co, m_max):
        
        #parámetros
        
        td = self.Td
        R = self.R
        fz = self.fz
        der_con = lambda c, km: 1 + km/c
        conc = np.zeros(len(R))
        
        for i in range(len(R)):
            
            ini = 0.000001
            b = -Co-km*math.log(Co) + R[i]*m_max[i]*fz[i]*td
            con = lambda c,b, km : c + km*math.log(c) + b
            it = 1
            while True:
                nv = ini - con(ini,b,km)/der_con(ini,km)
                if abs(ini - nv) < 0.0000001 or it > 50:
                    conc[i] = nv
                    Co = nv
                    break
                else:
                    ini = nv
                    it += 1
                    
        return conc
    
    
    def Cnit(self):
        
        self.C_NH4 = self.fconc( self.prop['Km_nit'], 60, self.krmax)
        self.C_NO3 = self.fconc( self.prop['Km_dnt'],1, self.Vmax)
        self.Ctotal = self.C_NO3 + self.C_NH4
        
        
            
    def graficar(self):
        x = self.prof
        y1 = self.C_NH4
        y2 = self.C_NO3
        y3 = self.Ctotal
        fig, ax = plt.subplots()
        ax.plot(x, y1)
        ax.plot(x, y2)
        ax.plot(x, y3)
        ax.legend(['CNH4','CNO3','CTotal'])
        plt.show()
        
    
    def newton(ini):
        
        con = lambda c : c + math.log(c) + 10 
        der_con = lambda c: 1 + 1/c
        it = 1
        while True:
            nv = ini - con(ini)/der_con(ini)
    
            if abs(ini - nv) < 0.0000001 or it>50:
                return nv
                break
            else:
                ini = nv
                it += 1
                

    
'''
################## Parámetros
D = 60.00	 # Soil depth 
dz = 0.80	 # Dz
################## 
HLR = 2.00#	HLR

aG = 0.03#	 aG
aVG = 0.02#	aVG
Ks= 14.75#	Ks
Yo = -8.61#	Yo
qr = 0.10#	qr
qs = 0.46#	qs 
n= 1.26#	n
m = 0.21#	m 
l = 0.50	#  l 
ho = 2.35	#  ho ##### Ajuste
cf = 1.00	#  cf #### Ajuste
#############
Co_NH4 = 60.00	#  Co NH4
Co_NO3 = 1.00	#  Co NO3

fs = 0.00#	  fs 
fwp = 0.00#	  fwp
swp = 0.15#	  swp
e2 = 2.27	 # e2
e3 = 1.10	#  e3 
kr_max = 56.00#	  kr max
Km_nit = 5.00#3	  Km,nit
bnit = 0.35#	  bnit 
sl = 0.67#	 sl
sh = 0.81#	  sh
##############
ednt = 3.77#	  ednt
Vmax = 2.56#	  Vmax 
Km_dnt = 5.00#	  Km,dnt
bdnt = 0.35#	  bdnt 
sdn = 0.00#	  sdn 
ac = 0.00#	  ac 
################
T = 18.50#	  T
kd = 1.46	#  kd
fr = 0.00	# fr
T_opt_nit = 25.00	 # Topt-nit
T_opt_dnt = 26.00	 # Topt-dnt
##############
ftdnt = 0.745861763178173	#  ftdnt Ajuste
ftnit = 0.687041502811397	 # ftnit Ajuste
rho = 1.50	  #rho
1.00	  #R'(NH4)
1.00	 # R (NO3)
MinCNH4 = 0.001	  #MinCNH4
MinCNO3 = 0.001	  #MinCNO3
'''
#Diccionario de propiedades
dic_h = {
'HLR':2.00,
'aG':0.025,#	 aG
'aVG':0.015,#	aVG
'Ks':14.75,#	Ks
'Yo':0,#	Yo
'qr':0.0980,#	qr
'qs':0.459,#	qs 
'n':1.26,#	n 
'l':0.50,	#  l 
######## nit
'swp':0.15,# swp ### sirve para calcular WFP
'fs':0.00,#	  fs 
'fwp':0.00,#	  fwp
'e2':2.27,	 # e2
'e3':1.10,	#  e3 
'kr_max':56.00,#	  kr max
'Km_nit':5.00,#3	  Km,nit
'bnit':0.35,#	  bnit 
'sl':0.67,#	 sl
'sh':0.81,#	  sh
####### dnt
'ednt':3.77,#	  ednt
'Vmax':2.56,#	  Vmax 
'Km_dnt':5.00,#	  Km,dnt
'bdnt':0.35,#	  bdnt 
'sdn':0.00,#	  sdn 
'ac':0.00,#	  ac 
##### R
'kd':1.46,
'rho':1.50	  #rho
}
##################### Modelo#################################
#perfil de profundidad
## creación del objeto modelo
mod = Modelo(60, 0.5, propiedades=dic_h)
#cálculo de profundidad
z = mod.z()
#Cálculo de presiones
y = mod.y = np.array(list(map(mod.psi,z)))
#perfil de humedad
hum = mod.tetha()
#saturación efectiva
se = mod.Se = np.array(list(map(mod.se,hum)))
#conductividad vertical
k = mod.K = np.array(list(map(mod.k,mod.Se)))
## función carbón
fz = mod.fcont_c()
#WFP
WFP = mod.fWFP()

## efecto de sw en la nit
fswnt = mod.fswnt()
## efecto de sw en la desn
fswdnt = mod.fswdnt()

### factor de retardación
R = mod.fR()

### kr max
krmax = mod.f_krmax()
##Vmax
Vmax = mod.fVmax()
## Concentraciones
mod.Cnit()
CO3 = mod.C_NO3
CNH4 = mod.C_NH4
##

mod.graficar()





