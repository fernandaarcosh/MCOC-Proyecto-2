# Proyecto 2 
from matplotlib.pylab import *

#Unidades base SI
_m = 1.
_mm = 1e-3*_m
_cm = 1e-2*_m

_kg = 1.
_gr = 1e-3*_kg
_ton = 1e3*_kg

_s = 1.

#velocidades iniciales
vfx = 10.0*_m/_s #m/s 
vfy = 0.1*_m/_s #m/s 

x0 = array([0.,1.*_mm], dtype=double)
v0 = array([1.,1.], dtype=double)
x0n = array([1.,1.*_mm], dtype=double)
v0n = array([1.,2.], dtype=double)

#velocidad y posicion actual
#xi = x0
#vi = v0

#velocidad y posicion en el instante mas 1
xim1 = zeros(2, dtype=double)
vim1 = zeros(2, dtype=double)

g = 9.81 *_m/_s**2
d= 1*_mm
rho = 2650.*_kg/(_m**3)
Cd = 0.47 # particula esferica
rho_agua = 1000*_kg/(_m**3)

V = (4./3.)*pi*(d/2)**3
A = pi*(d/2)**2
m = rho*V                 # masa de la particula

#Inicializar Euler en x0

dt= 0.001*_s   #paso de tiempo
tmax = 2*_s   # tiempo maximo de simulacion


ti = 0.0*_s   #tiempo actual

W = array([0, -m*g])
fB = array([0., rho_agua*V*g])
vf = array([vfx,vfy])

norm = lambda v: sqrt(dot(v,v))

t = arange(0., tmax, dt)
Nt = int32(2*tmax /dt) 

k_penal = 1000*0.5*Cd*rho_agua*A*norm(v0)/(1*_mm)

#x_store = zeros((2,Nt))
#v_store = zeros((2,Nt))
#t_store = zeros((Nt))

# Empieza metodo de euler

def particula(z,t):
	k = 0

 	zp = zeros(len(z))
 	for a in range((len(z)/4)):
 		xi = z[k:2+k]
 		vi = z[2+k:4+k]
 		vf = array([vfx,vfy])
 		vrel = vf - vi
 		fD = (0.5*Cd*rho_agua*norm(vrel)*A)*vrel  # Fuerzza de Drag
 		Fi = W + fD + fB          # Sumatoria de todas las fuerzas

	if xi[1] < 0:
		Fi[1] += -k_penal*xi[1]
	
		zp[k:2+k] = vi

		zp[2+k:4+k] = Fi/m

	k += 4
	return zp

from scipy.integrate import odeint

# Para un numero deseado de particulas
#n = input("Numero de Particulas") # numero de particulas
#x_store = []
#v_store = []
#for i in range (0,n):
#	p1 = float(random.randit(0,n))
#	x0n = array([0., p1*_mm], dtype=double)
#	v0n = array([1., 1.], dtype=double)
#	xin = x0n
#	vin = v0n
# 	xin1 = zeros(2, dtype=double)
# 	vin1 = zeros(2, dtype=double)
#	z0 = zeros(4)
#	z0[:2] = x0n
#	z0[2:] = v0n
#	z = odeint(particula, z0, t)
#	xn = z[:, :2]
#	vn = z[:, 2:]

#x_store.append(xn)
#v_store.append(vn)


	
#Para 2 particulas
n = 2
z0 = zeros(4*n)
z0[:2] = x0
z0[2:4] = v0
z0[4:6] = x0n
z0[6:8] = v0n

z = odeint(particula, z0, t)

x = z[:,:2]
v = z[:,2:4]
xn = z[:,4:6]
vn = z[:,6:8]

figure()
plot(x[:, 0.], x[:, 1.], label="x")
plot(xn[:,0.], xn[:, 1.], label="xn")
ylim([0, 10*_mm])
plt.legend()

figure()
subplot(2, 1, 1)
plot(t, x[:, 0.], label="x")
plot(t, x[:, 1], label="y")
plt.legend()

subplot(2,1,2)
plot(t, v[:, 0.], label="vx")
plot(t, v[:, 1], label="vy")

figure()
subplot(2, 1, 1)
plot(t, xn[:, 0.], label="xn")
plot(t, xn[:, 1], label="yn")
plt.legend()

subplot(2,1,2)
plot(t, vn[:, 0.], label="vxn")
plot(t, vn[:, 1], label="vyn")

show()
