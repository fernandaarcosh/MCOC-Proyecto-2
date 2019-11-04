from matplotlib.pylab import *
from parameters import *

norm = lambda v: sqrt(dot(v,v))

def propiedades_area_volumen_masa(d):
	area= pi * (d/2) ** 2
	vol = (4./3.) * pi * (d/2) ** 3
	masa = rho_particula * vol
	return area, vol, masa

def fuerzas_hidrodinamicas(x,v,d,area,masa):

	xtop = x + (d/2)*jhat
	xbot = x - (d/2)*jhat
	vf = velocity_field(x + 0*jhat)

	vrelf_top = abs(velocity_field(xtop)[0])
	vrelf_bot = abs(velocity_field(xbot)[0])

	vrel = vf - v

	Cd = 0.47
	fD = (0.5*Cd*alpha*rho_agua*norm(vrel)*area)*vrel

	fL = (0.5*CL*alpha*rho_agua*(vrelf_top - vrelf_bot)*area)*vrel[0]*jhat
	fW = (-masa*g)*jhat

	Fh = fW + fD + fL

	return Fh

def fuerza_impacto_suelo(x,v,d):
	N = around(x[0]/d)
	r = x - (N * d) * ihat
	delta = norm(r)- d
	if delta < 0:
		n = r/ norm(r)
		Fi = -k_penal * delta * n
	else:
		Fi = 0. * r 
	return Fi






def zp_una_particula(z,t,d=d):
	zp = zeros(4)

	x = z[0:2]
	v = z[2:4]

	area, vol, masa = propiedades_area_volumen_masa(d)
	Fh = fuerzas_hidrodinamicas(x,v,d,area,masa)
	Fi = fuerza_impacto_suelo(x,v,d)

	sumaF = Fh + Fi

	zp[0:2] = v
	zp[2:4] = sumaF/masa

	return zp



def zp_todas_las_particulas(z,t):
	zp = zeros(4 * Nparticulas)

	for i in range(Nparticulas):
		di = d
		xi = z[4 * i: (4 * i + 2)]
		vi = z[4 * i + 2 : (4*i +4)]

		area,vol,masa = propiedades_area_volumen_masa(d)
		Fh = fuerzas_hidrodinamicas(x,v,d,area,masa)
		Fi = fuerza_impacto_suelo(x,v,d)

		sumaF = Fh + Fi

		zp[4*i:(4*i+2)] = vi
		zp[4*i+2:(4*i+4)] = sumaF/masa   # decia m

	zp += zp_choque_M_particula(z,t,M = Nparticulas)

	return zp


def zp_M_particulas(z,t,M):
	zp = zeros(4*M)

	for i in range(M):
		di = d
		zi = z[4*i:(4*i+4)]
		vi = z[4*i+2:(4*i+4)]

		zp[4*i:(4*i+4)] = zp_una_particula(zi,t,di)

	zp += zp_choque_M_particula(z,t,M=M)

	return z|p


def zp_choque_M_particulas(z,t,M):
	zp = zeros(4*M)
	for i in range(M):
		xi = z[4*i:(4*i+2)]
		di = d
		area_i,vol_i,masa_i = propiedades_area_vol_masa(di)
		for j in range(i+1,M):
			xj = z[4*j:(4*j+2)]
			dj = d
			rij = xj - xi
			norm_rij = norm(rij)
			if norm_rij < 0.5*(di+dj):
				area_j,vol_j,masa_j = propiedades_area_vol_masa(dj)
				delta = 0.5 * (di+dj) - norm_rij
				nij = rij / norm_rij
				Fj = k_penal * delta * nij
				Fi = -Fj
				zp[4*i+2:(4*i-4)] += Fi/masa_i
				zp[4*j+2:(4*j+4)] += Fj/masa_j
	return zp
