from matplotlib.pylab import *

norm = lambda v: sqrt(dot(v,v))
#Completo

Nparticulas = 2


#Unidades base son SI
_m = 1.
_kg = 1.
_s = 1.
_mm = 1e-3 *_m
_cm = 1e-2*_m
_gr = 1e-3 * _kg
_in = 2.54*_cm

g = 9.81*_m/_s**2
d = 0.15*_mm

rho_agua = 1000.*_kg/(_m**3)
rho_particula = 2650.*_kg/(_m**3)

dt = 0.00001*_s    	# paso de tiempo 
tmax = 0.05*_s  # tiempo maximo de simulacion

Cd = 0.47 # Para particula esferica
Cm = 0.5
CL = 0.2
Rp = 73.

ustar = 0.14


tau_star = 0.067  


R = (rho_particula/rho_agua - 1)
alpha = 1/(1 + R + Cm)

ihat = array([1,0])
jhat = array([0,1])

tau_cr = 0.22*Rp**(-0.6)+0.06*10**(-7*Rp**(-0.6))   # tau critico
ustar = sqrt(tau_star * g * Rp * d)

print "tau_star=",tau_star
print "tau_cr=",tau_cr
print "tau_star/tau_cr=",tau_star/tau_cr
print "ustar=",ustar




def velocity_field(x):
	z = x[1] /d

	if z > 1./30:
		uf = ustar*log(30.*z)/0.41
		uf = uf * (uf>0)
	else:
		uf = 0
		
	return array([uf,0])
	

vfx = velocity_field([0, 10*d])[0]
A = pi*(d/2)**2
k_penal = 0.5*Cd*rho_agua*A*norm(vfx)**2/(d/20)
