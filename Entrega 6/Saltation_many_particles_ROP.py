from matplotlib.pylab import *
import time 

norm = lambda v: sqrt(dot(v, v)) #norma vector
Nparticulas = 2

#Unidades SI
_m = 1.
_kg = 1.
_s = 1.
_mm = 1e-3 * _m
_gr = 1e-3 * _kg

g = 9.81*_m/_s**2 #gravedad
d = 0.15 * _mm #diametro particula 

rho_agua = 1000.*_kg/(_m**3) #densidad agua
rho_particula = 2650.*_kg/(_m**3) #densidad arena

dt = 0.00001*_s #intervalo tiempo
tmax = 0.05*_s #tiempo maximo

Cd = 0.47 #coeficiente de drag
Cm = 0.5 
CL = 0.2 #coeficiente de lif
Rp = 73.

tau_star = 0.067

R = (rho_particula/rho_agua -1)
alpha = 1/(1 + R + Cm)
ihat = array([1,0])
jhat = array([0,1])

tau_cr = 0.22*Rp**(-0.6)+0.06*10**(-7*Rp**(-0.6))   #tau critico
ustar = sqrt(tau_star * g * Rp * d)

print "tau_star=", tau_star
print "tau_cr=", tau_cr
print "tau_star/tau_cr=", tau_star/tau_cr
print "ustar=", ustar

def velocity_field(x):
	z = x[1] / d

	if z > 1./30:
         uf = ustar*log(30.*z)/0.41
         uf = uf * (uf>0)
	else:
		uf = 0
		
	return array([uf,0])
	
vfx=velocity_field([0, 10*d])[0]
A=pi*(d/2)**2
k_penal=0.5*Cd*rho_agua*A*norm(vfx)**2/(d/20)

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

	# delta = x[1]-d/2
	# if delta < 0:
	#	Fi += -k_penal*delta*jhat
	# return Fi


def zp_una_particula(z,t,d=d):  #sumatoria de fuerzas para una particula para cuando no hay choque entre particulas
	zp = zeros(4)
	x = z[0:2]
	v = z[2:4]
	area, vol, masa = propiedades_area_volumen_masa(d)
	Fh = fuerzas_hidrodinamicas(x,v,d,area,masa)
	Fi = fuerza_impacto_suelo(x,v,d)
	sumaF = Fh + Fi
	zp[0:2] = v
	zp[2:4] = sumaF/masa
	return zp  #nueva velocidad y aceleracion



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
		zp[4*i+2:(4*i+4)] = sumaF/masa   
	zp += zp_choque_M_particula(z,t,M = Nparticulas)
	return zp

def zp_choque_M_particulas(z,t,M): #funcion de sumatorias de fuerzas, se le suma a cada particula su choque con las demas particulas
	zp = zeros(4*M)
	for i in range(M):
		xi = z[4*i:(4*i+2)]
		di = d
		area_i,vol_i,masa_i = propiedades_area_volumen_masa(di)
		for j in range(i+1,M):
			xj = z[4*j:(4*j+2)]
			dj = d
			rij = xj - xi
			norm_rij = norm(rij)
			if norm_rij < 0.5*(di+dj):
				area_j,vol_j,masa_j = propiedades_area_volumen_masa(dj)
				delta = 0.5 * (di+dj) - norm_rij
				nij = rij / norm_rij
				Fj = k_penal * delta * nij
				Fi = -Fj
				zp[4*i+2:(4*i-4)] += Fi/masa_i #aceleracion de cada particula
				zp[4*j+2:(4*j+4)] += Fj/masa_j
	return zp


def zp_M_particulas(z,t,M): #para ver la interaccion entre particulas
	zp = zeros(4*M)
	for i in range(M):
		di = d
		zi = z[4*i:(4*i+4)]
		vi = z[4*i+2:(4*i+4)]
		zp[4*i:(4*i+4)] = zp_una_particula(zi,t,di)
	zp += zp_choque_M_particulas(z,t,M=M)
	return z|p

#reuse_initial_condition = True
reuse_initial_condition = False

doit = True
#doit = False

start = time.time ()

tiempo_bloque_1 = 0
tiempo_bloque_2 = 0
t = arange (0,tmax,dt)
Nt = len(t)

Nparticulas = 80

if reuse_initial_condition: 
	print "Reusing initial conditions"
	data = load("initial_condition.npz")
	x0 = data ["x0"]
	y0 = data ["y0"]
	vx0 = data ["vx0"]
	vy0 = data ["vy0"]
	Nparticulas = data ["Nparticulas"]
else: 
	print "Generating new initial conditions"	
	itry = 1
	while True:
		dmin = infty
		x0 = 800*d*rand(Nparticulas)
		y0 = 5*d*rand(Nparticulas) + 1*d
		for i in range (Nparticulas):
			xi, yi = x0[i], y0[i]
			for j in range (i+1, Nparticulas):
				xj, yj = x0[j], y0[j]
				dij = sqrt ((xi - xj)**2 + (yi - yj)**2)
				dmin = min (dmin, dij)
		print "Try # " , itry, "dmin/d =", dmin/d	
		if dmin > 0.9*d:
			break
		itry += 1

	vx0 = ustar*rand(Nparticulas)		
	vy0 = 0
	savez ("initial_condition.npz", x0=x0, y0=y0, vx0=vx0, vy0=vy0, Nparticulas=Nparticulas)

t = arange (0, tmax, dt)
Nt = len (t)

from scipy.integrate import odeint 

# z = zeros ((Nt, 4*Nparticulas))
zk = zeros ((4*Nparticulas))
zkm1 = zeros ((4*Nparticulas))

zk[0::4] = x0 
zk[1::4] = y0
zk[2::4] = vx0
zk[3::4] = vy0



#print Nt
#exit (0)
fout = open("resultado.txt", "w")

done = zeros (Nparticulas,dtype=int32)
impacting_set = zeros(Nparticulas, dtype=int32)

print "Integrando"
k = 0

if doit:
	while dt*k < int(tmax/dt-1)*dt: 

		
		fout.write("{} ".format(dt*k)) #escribo el paso actual de tiempo


		ti = time.time()				
		savetxt(fout, zk, fmt='%.5e ', newline=" ") 
		#fout_z[k,0] = dt * k
		#fout_z[k,1:] = zk
		tf = time.time()

		tiempo_bloque_1 += tf-ti

		fout.write("\n")

		if k % 100 == 0:
			print "k = {}    t = {}  ".format(k, k*dt)
		done *= 0

		ti = time.time()

			
		for i in range (Nparticulas):
			irange = slice(4*i, 4*i+4)
			
			zk_i = zk[irange]
			

			di = d
			if done[i] == 0:

				hay_impacto = False
				impacting_set *=0
				M = 1

				for j in range(i+1, Nparticulas):
					jrange = slice(4*j, 4*j+4)
					zk_j = zk[jrange]

					dj = d
					rij = zk_j[0:2] - zk_i[0:2]

					if norm (rij) < 0.5*(di-dj)*3:
						hay_impacto = True
						impacting_set[0] = i
						impacting_set[M] = j
						M+=1

				if hay_impacto:
				
					zk_all = zk_i
					for j in impacting_set[1:M]:
						jrange = slice(4*j, 4*j+4)
						zk_j = zk[jrange]
						zk_all = hstack((zk_all, zk_j))

					ti = time.time()
					
					zkm1_all = odeint(zp_M_particulas, zk_all, [dt*k, dt*(k+1)], args=(M,))

					#tf = time.time()
					#tiempo_bloque_1 += tf - ti

					zkm1 [irange] = zkm1_all[1,0:4]
 
					done[i] = 1
					pos_j = 1
					for j in impacting_set[1:M]:
						jrange = slice (4*j, 4*j+4)
						zkm1[jrange] = zkm1_all[1,4*pos_j:4*pos_j+4]
						done[j] = 1
						pos_j += 1		
				
				else:
					ti = time.time()

					zkm1_i = odeint(zp_una_particula, zk_i, [dt*k, dt*(k+1)])

					tf = time.time()
					tiempo_bloque_2 += tf - ti

					zkm1[irange] = zkm1_i [1,0:4]
					done [i] = 1

		tf = time.time()

		tiempo_bloque_2 += tf-ti 
					
		zk = zkm1
		k += 1			

end = time.time()

print "tiempo bloque 1: ", tiempo_bloque_1
print "tiempo bloque 2: ", tiempo_bloque_2
print "tiempo total: ", end - start		
print "Tiempo de escritura es {}% del tiempo total.". format (tiempo_bloque_1/(end-start)*100)			

fout.close()