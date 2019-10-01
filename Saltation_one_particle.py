# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 11:38:57 2019

@author: aniba
"""

from numpy import *
from matplotlib.pylab import *

#unidades base
_m = 1
_kg = 1
_s = 1
_mm = 1e-3*_m
_gr = 1e-3*_kg

vfx = 5.0*_m/_s#m/s
vfy = 0.0*_m/_s

x0 = array([0., 1.], dtype=double)
v0 = array([1., 1.], dtype=double)

xi = 0 #posicion actual
vi = 0 #velocidad actual
xim1 = zeros(2, dtype=double)#posicion siguiente
vim1 = zeros(2, dtype=double)#velocidad siguiente

g = 9.81*_m/(_s**2)
d = 1.*_mm
rho = 2700*_kg/(_m**3)
m = rho*(4./3./8.)*pi*(d**3) #masa de la particula
Cd= 0.47 #coef de drag para un apaticula

#Inicialzar Euler en x0

dt = 2e-6*_s #s, paso de tiempo
tmax = 1*_s #tiempo max de simulacion
ti = 0*_s #tiempo actual

W = array([0, -m*g])
vf = array([vfx, vfy])

Nt = int32(2*tmax / dt)
x_store = zeros((2,Nt))
v_store = zeros((2,Nt))
t_store = zeros(Nt)

#Metodo de Euler
i = 0
while ti < tmax:
    
    #evaluar v. realtiva
    if i % 100 == 0:
        print("ti = ", ti  ,"  |xi| = ", sqrt(dot(xi, xi)))
    vrel = vf - vi
    norm_vrel = sqrt(dot(vrel,vrel))
    
    #evaluar fzas sobre la particula
    fD = 0,5*Cd*norm_vrel*vrel
    Fi = W + fD
    
    #evualuar aceleracion
    ai = Fi / m
    
    #integrar
    xim1 = xi + vi * dt + ai*(dt**2/2)
    vim1 = vi + ai*dt
    
    #avanzar a sgte paso
    x_store[:, i] = xi
    v_store[:, i] = vi
    t_store[i] = ti
    
    ti += dt
    i += 1
    xi = xim1
    vi = vim1


#guardar tu ultimo paso
x_store[:,i] = xi
v_store[:,i] = vi
t_store[i] = ti

print(x_store)

figure()
plot(x_store[0,:i], x_store[1, :i])
show()

