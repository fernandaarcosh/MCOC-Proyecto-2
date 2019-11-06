from matplotlib.pylab import * 
from parametros import *
import time
import h5py
import scipy as sp 

norm = lambda v: sqrt(dot(v,v))

# Funcion que determina las propiedades fisicas de la particula:
def propiedades_area_volumen_masa(d):
	area= pi * (d/2) ** 2
	vol = (4./3.) * pi * (d/2) ** 3
	masa = rho_particula * vol
	return area, vol, masa

# Funcion fuerzas hidrodinamicas
def fuerzas_hidrodinamicas(x,v,d,area,masa):

	xtop = x + (d/2)*jhat
	xbot = x - (d/2)*jhat
	vf = velocity_field(x + 0*jhat)

	vrelf_top = abs(velocity_field(xtop)[0])
	vrelf_bot = abs(velocity_field(xbot)[0])

	vrel = vf - v

	Cd = 0.47
	fD = (0.5*Cd*alpha*rho_agua*norm(vrel)*area)*vrel # Fuerza de Drag

	fL = (0.5*CL*alpha*rho_agua*(vrelf_top - vrelf_bot)*area)*vrel[0]*jhat # Fuerza de Lift
	
	fW = (-masa*g)*jhat # Fuerza de Peso

	Fh = fW + fD + fL
	return Fh

# Funcion que determina la fuerza de impacto que se tiene en el suelo
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

# Funcion que configura los datos de posicion, velocidad y fuerzas sobre y de una particula
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


def M_particulas(z,t,M):
	zp = zeros(4*M)

	for i in range(M):
		di = d
		zi = z[4*i:(4*i+4)]
		vi = z[4*i+2:(4*i+4)]

		zp[4*i:(4*i+4)] = zp_una_particula(zi,t,di)

	zp += zp_choque_M_particula(z,t,M=M)

	return zp

# Funcion que determina el choque entre particulas, es decir cuales chocan y cuales no.
def choque_M_particulas(z,t,M):
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


#reuse_inicial_condicion = True
reuse_initial_condition = False

doit = True
#doit = False
 
Tiempo_Inicio= time.time()
tiempo_bloque_1=0
tiempo_bloque_2=0
t=arange(0,tmax,dt)
Nt=len(t)

# Archivo de texto que guardara los datos de la simulacion realizada
fout=open("Resultados.txt","w")
#fout = open("resultado.txt","w")

if reuse_initial_condition:
	print "reusing_initial_conditions"
	data = load("initial_condition.npz")
	x0 = data["x0"]
	y0 = data["y0"]
	vx0 = data["vx0"]
	vy0 = data["vy0"]
	Nparticulas = data["Nparticulas"]
else:
	print "Generando nuevas condiciones iniciales"
	itry = 1
	while True:
		dmin = infty																					
		x0 = 800 * d *rand(Nparticulas)
		y0 = 5 * d *rand(Nparticulas) + 1 * d								
		for i in range(Nparticulas):
			xi,yi = x0[i], y0[i]
			for j in range(i+1, Nparticulas):
				xj,yj = x0[j],y0[j]
				dij = sqrt((xi-xj)**2 +(yi-yj)**2)
				dmin = min(dmin,dij)												
		print "try #", itry, "dmin/d = ",dmin/d	
		if dmin >  0.9 * d:
			break
		itry += 1

	vx0 = ustar * rand(Nparticulas)
	vy0 = 0 
	savez("initial_conditions.npz",x0=x0,y0=y0,vx0=vx0,vy0=vy0,Nparticulas=Nparticulas)


t = arange(0,tmax,dt)
Nt = len(t)



from scipy.integrate import odeint


zk = sp.zeros((4*Nparticulas))
zkm1 = sp.zeros((4*Nparticulas))

zk[0::4] = x0 
zk[1::4] = y0
zk[2::4] = vx0
zk[3::4] = vy0

# Se pasan los datos a binario para ahorrar memoria y optimizar el codigo
fout = h5py.File("binario", "w")
fout_parametros = fout.create_group("Parametros")
fout_parametros["dt"] = dt
fout_parametros["g"] = g
fout_parametros["d"] = d
fout_parametros["rho_agua"] = rho_agua
fout_parametros["rho_particula"] = rho_particula
fout_parametros["tmax"] = tmax
fout_parametros["Cd"] = Cd
fout_parametros["Cm"] = Cm
fout_parametros["CL"] = CL
fout_parametros["Rp"] = Rp
fout_parametros["ustar"] = ustar
fout_parametros["tau_star"] = tau_star
fout_parametros["R"] = R
fout_parametros["alpha"] = alpha
fout_parametros["ihat"] = ihat
fout_parametros["jhat"] = jhat
fout_parametros["tau_cr"] = tau_cr
fout_parametros["A"] = A
fout_parametros["k_penal"] = k_penal
fout_z = fout.create_dataset("z",(Nt, 1+ 4*Nparticulas), dtype=double)

done = zeros(Nparticulas, dtype = int32)
impacting_set = zeros(Nparticulas,dtype = int32)

print "Integrando"
k = 0
if doit:
	while dt*k < int(tmax/dt-1)*dt:  

		#guarda el tiempo y los suma para el bloque de particulas que no chocan
		tiempo_actual = time.time()
		fout_z[k,0] = dt*k
		fout_z[k,1:] = zk
		tiempo_iteracion = time.time()
		tiempo_bloque_1 += tiempo_iteracion - tiempo_actual


		if k % 100 == 0:
			print "k = {} t = {} ".format(k,k*dt)

		done *= 0


		#guarda el tiempo y los suma para el bloque de particulas que chocan
		tiempo_actual = time.time()

		for i in range(Nparticulas):
			irange = slice(4*i, 4*i+4)

			zk_i = zk[irange]
			di = d
			if done[i] == 0:
				hay_impacto = False
				impacting_set *= 0
				M = 1
				for j in range(i+1,Nparticulas):
					jrange = slice(4*j, 4*j+4)														
					zk_j =zk[jrange] 													
					dj = d
					rij = zk_j[0:2] - zk_i[0:2]
					if norm(rij) < 0.5 * (di-dj) * 3: #problemas con norm
						hay_impacto = True
						impacting_set[0] = i
						impacting_set[M] = j 
						M += 1

				if hay_impacto:

					zk_all = zk_i
					for j in impacting_set[1:M]:											
						jrange = slice(4*j , 4*j+4)										
						zk_j = zk[jrange]
						zk_all = hstack((zk_all, zk_j))										

					zkm1_all =odeint(M_particulas,zk_all, [dt*k, dt*(k+1)],args=(M,))
			
					zkm1[irange] = zkm1_all[1,0:4]

					done[i] = 1
					pos_j = 1
					for j in impacting_set[1:M]:
						jrange = slice(4*j, 4*j + 4)
						zkm1[jrange] = zkm1_all[1,4*pos_j:4*pos_j+4]
						done[j] = 1
						pos_j +=1

				else:
					
					zkm1_i = odeint(zp_una_particula, zk_i, [dt*k, dt*(k+1)]) ###cambien - por +								

					zkm1[irange] = zkm1_i[1,0:4]
					done[i] = 1

		tiempo_iteracion = time.time()
		tiempo_bloque_2 += tiempo_iteracion - tiempo_actual

		zk=zkm1
		k+=1

Tiempo_Final = time.time()
fout.close()

print 'Tiempo bloque 1, particulas que no chocan: ',tiempo_bloque_1
print 'Tiempo bloque 2, particulas que chocan: ',tiempo_bloque_2
print 'Tiempo Total:', Tiempo_Final - Tiempo_Inicio

with h5py.File("binario", 'r') as f:
    print("Keys: %s" % f.keys())
    a_group_key = list(f.keys())[1]
    data = list(f["z"])

d = 0.15e-3 

Nparticulas = (len(data[0]) -1) /4

figure()

color = "006B93"
plt = gca()
colorlist = []
xi =[]
yi =[]

#para cada particula
for i in range(Nparticulas):
	#para cada particula guarda la trayectoria
	for j in range(len(data)-1):
		xi.append(data[j][1 + 4*i] / d)
		yi.append(data[j][1 + 4*i + 1] / d)
	
	col=rand(4)
	colorlist.append(col)
	plt.plot(xi[0::100],yi[0::100],"o",color=col)
	plt.plot(xi,yi,"--",color=col,alpha=0.6)
	xi = []
	yi = []

plt.set_ylim([0,8])
plt.axhline(0.,color="b",linestyle="--")
plt.axhline(1/30.,color="gray",linestyle="--")
plt.set_xlabel("${x}/{d}$") # Unidades eje x
plt.set_ylabel("${z}/{d}$")	# Unidades eje y
if Nparticulas == 1:
	titulo_grafico = 'Trayectoria de {} particulas'.format(Nparticulas)
else:
	titulo_grafico = 'Trayectoria de {} particulas'.format(Nparticulas)
plt.set_title(titulo_grafico)
tight_layout()

nombre_imagen = 'Grafico {} particulas'.format(Nparticulas)
savefig(nombre_imagen)

show()
