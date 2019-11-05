from matplotlib.pylab import *
from parameters import *
from funciones import *
from time import time


reuse_initial_condition = True
# reuse_initial_condition = False

doit = True
# doit = False

inicio = time.time()

tiempo_bloque_1 = 0
tiempo_bloque_2 = 0
t = arrange(0, tmax, dt)
Nt = len(t)

Nparticulas = 2

if reuse_initial_condition:
	print "Reusando condiciones iniciales"
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
		x0 = 670*d*rand(Nparticulas)
		y0 = 5*d*rand(Nparticulas) + 1*d
		for i in range(Nparticulas):
			xi, yi = x0[i], y[i]
			for j in range(i + 1, Nparticulas):
				xj, yj = x0[j], y0[j]
				dij = sqrt((xi - xj)**2 - (yi -yj)**2)
				dmin = min(dmin, di)
		print "Try #", itry, "dmin/d= ", dmin/d
		if dmin > 0.9*d:
			break
		itry += 1

	vx0 = ustar*rand(Nparticulas)
	vy0 = 0
	savez("initial_condition.npz", x0 = x0, y0 = y0, vx0 = vx0, vy0 = vy0, Nparticulas = Nparticulas)


def particula(z,t):
		zp = zeros(4*Nparticulas)

		for i in range(Nparticulas):
			di = d 
			xi = z[4*i:(4*i+2)]
			vi = z[4*i+2:(4*i+4)]

			xtop = xi + (d/2)*jhat
			xbot = xi - (d/2)*jhat
			vf = velocity_field(xi + 0*jhat)
			vf_top = abs(velocity_field(xtop)[0]) 
			vf_bot = norm(velocity_field(xbot)[0]) 
			vrel = vf - vi
			fD = (0.5*Cd*alpha*rho_agua*norm(vrel)*A)*vrel  #formula wiki
			#fD = alpha*(R*(d*g/(ustar**2))-(3./4.)*Cd*(vrel)*norm(vrel)) # formula PM
			fL = (0.5*CL*alpha*rho_agua*norm(vf_top**2 - vf_bot**2)*A)*jhat #formula wiki
			#fL = alpha*(3/4*CL*(norm(vf_top)**2 - norm(vf_bot)**2)) # formula PM
			fB =  alpha*(-rho_agua*g*V*jhat) # fromula de empuje


			Fi = W + fD + fL + fB

			xs = (xi[0]//d)*d + d/2
			cxs = array([xs,0])
			dxs = xi - cxs
			if norm(dxs) < d:
				delta = d - norm(dxs)
				nxs = dxs/norm(dxs)
				Fi += k_penal*delta*nxs

			zp[4*i:(4*i+2)] = vi
			zp[4*i+2:(4*i+4)] = Fi/m

			for i in range(Nparticulas):
				xi = z[4*i:(4*i+2)]
				for j in range(Nparticulas):
					if i > j:
						xj = z[4*j:(4*j+2)]
						rij = xj -xi
						if norm(rij) < d:
							delta = d - norm(rij)
							nij = rij/norm(rij)
							Fj = k_penal*delta*nij
							Fi = -k_penal*delta*nij
							zp[4*i+2:(4*i+4)] += Fi/m
							zp[4*j+2:(4*j+4)] += Fj/m

		return zp

#Choque entre particulas:
# j = -f*(1+e)*n*(vxi-vxj)*(m/2)

#choque con el suelo
# r/d = 0.5*(cos(tethab)-tan(tethain)*sin(tethab))
#promedio de los saltos de las particulas, H/d para un Tstar/tcritico definido.
# caracteristicas de pc, raficar tiempo de particulas respecto a cuantose demora.

from scipy.integrate import odeint
zk= zeros((4*Nparticulas))
zkm1 = zeros((4*Nparticulas))

zk[0::4] = x0
zk[1::4] = y0
zk[2::4] = vx0
zk[3::4] = vy0


arch_part = open("resultado.txt", "w")


#z=odeint(particula,z0,t)
hecho = zeros(Nparticulas, dtype=int32)
impacto = zeros(Nparticulas, dtype=int32)

print "Integrando"
k = 0

if yahecho:
	while dt*k < int(tmax/dt-1)*dt:
		
		arch_part.write("{}".format(dt*k))
		savetxt(arch_part, zk, fmt= "%.18e ", newline = "")
		arch_part.write("\n")

		if k % 100 == 0:
			print "k = {}  t = {}" .format(k, k*dt)
		#zk = z[k, :]
		hecho *= 0
		for i in range(Nparticulas):
			irange = slice(4*i, 4*i + 4)
			zk_i = zk[irange]
			di = d
			if hecho[i] == 0:
				hay_impacto = False
				impacto *= 0
				M = 1
				for j in range(i+1, Nparticulas):
					jrange = slice(4*j, 4*j +4)
					zk_j = zk[jrange]
					dj = d
					rij = zk_j[0:2]- zk_i[0:2]
					if norm(rij) < 0.5*(di + dj)*3:
						hay_impacto = True
						impacto[0] = i
						impacto[M] = j
						M += 1
				if hay_impacto:
					zk_all = zk_i
					for j in impacto[1:M]:
						xjrange = slice(4*j, 4*j+4)
						zk_j = zk[jrange]
						zk_all = hstack((zk_all, zk_j))
					ti = time.time()

					zkm1_all = odeint(zp_M_particulas,zk_all)

					tf = time.time()
					tiempo_bloque_1 += tf - ti

					z[k+1, irange] = zkm1_all[1, 0:4]
					hecho [i] = 1
					post_j = 1
					for j in impacto[1:M]:
						jrange = slice(4*j, 4*j+4)
						z[k+1, jrange] = zkm1_all[1, 4*post_j:4]
						hecho[j] = 1
						post_j += 1
				else:
					ti = time.time()

					zkm1_i = odeint(zp_una_particula, zk_i, [dt*k, dt*(k+1)])

					tf = time.time()
					tiempo_bloque_2 += tf - ti

					# z[k+1, irange] - zkm1_i[1, 0:4]
					zkm1 = zkm1_i[1, 0:4]
					hecho[i] = 1
		k += 1

#else:
#	data = load("solution.npz")
#	t = data["t"]
#	z = data["z"]
#	dt = data["dt"]

fin = time.time()

print "Tiempo Bloque 1: ", tiempo_bloque_1
print "Tiempo Bloque 2: ", tiempo_bloque_2
print "Tiempo Total: ", fin - inicio

fig=figure()

ax=gca()
for i in range(Nparticulas):
	xi=z[:,4*i]/d
	yi=z[:,4*i+1]/d
	col= rand(3)
	plot(xi[0],yi[0],"o",color="r")

	plot(xi,yi,color=col)

xmax = 0
for a in range(Nparticulas):
	xi=z[:,4*i:]
	xmax1 = max(xi[:,0])
	xmax = max([xmax,xmax1])

l = (xmax/d)*100
x = linspace(0, xmax, l)
x_mod_d = (x % d) - d/2
y = sqrt((d/2)**2 - x_mod_d**2)	
plot(x/d, y/d)

ax.axhline(d/2,color="k", linestyle="--")
ax.axhline(0,color="k", linestyle="--")
plt.xlabel("Avance direccion X (mm)")
plt.ylabel("Altura direccion Y (mm)")
plt.title("Movimiento de particulas (plano XY)")
plt.legend()

tiempo_final= time() 

tiempo_de_compilacion= tiempo_final-tiempo_inicial	
print 'El tiempo de ejecucion fue:',tiempo_de_compilacion, "segundos" #En segundos

#axis("equal")

show()

