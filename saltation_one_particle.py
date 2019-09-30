
# Para simular el salto de una particula, utilizamos la Free flight particle equation, ignorando la fuerza
# de Magnus y la Fuerza historica o de Basset

#------------------------------------------
#Formula perfil de velocidad

# Uf(x)= 1/0.41 * ln(30*x)

# A la funcion uno le entrega la ubicacion en el eje x

#--------------------------------------------

#Datos

#volumen_particula= 8.80*(10**-9) # [m3]
#diametro_particula= 2.1 # [mm]
#masa_particula= 2.33*(10**-2) # [gr]

#radio_particula= diametro_particula/2

#phro= masa_particula/volumen_particula	

#coeficiente_drag= 0.47 # se asume particula esferica

#coeficiente_lift= (2*radio_particula)/(phro*(volumen_particula**2)*A)

#g= 9.8 #[m/s2] ---> gravedad
#-------------------------------------------------------------------------------------------------

# Proyecto 2 MCOC 
from matplotlib.pylab import *

# Unidades base SI (m, kg, s)
_mm = 1e-3*_m
_m = 1.
_km = 1e3*_m

_gr = 1e-3*_kg 
_kg = 1.
_ton = 1e3*_kg

_s = 1.
_min = 60*_s
_hr = 3600*_s

# Variables 
vf_x = 5.0*_m/_s
vf_y = 0.0*_m/_s

x_0 = array([0., 1.], dtype=doube)
v_0 = array([1., 1.], dtype=double)

x_i = zeros(2, dtype=double)      # posicion actual
v_i = zeros(2, dtype=double)      # velocidad actual
x_i_m1 = zeros(2, dtype=double)   # posicion siguiente
v_i_m1 = zeros(2, dtype=double)   # velocidad siguiente

g = 9.81*_m/(_s**2)         # gravedad
d = 2.1*_mm                   # diametro de la particula
rho = 1550*_kg/(_m**3)      # densidad de la particula, considerando que sea arena 
m = rho*pi*((d/2)**3)*(4/3) # masa de la particula
Cd = 0.47                   # coeficiente de Drag para particula esferica

# Inicializar Euler en x_0
dt = 2e-6*_s    # paso de tiempo 
t_max = 1*_s    # tiempo máximo de simulacion
t_i = 0.*_s     # tiempo actual

# Fuerzas
W = array([0., -m*g])  #Fuerza de peso 

# velocidades
v_f = array([vf_x, vf_y])

# Para ir guardando
Nt = int32(t_max/dt)
x_store = x_0 #zeros((2, Nt))
v_store = v_0 #zeros((2, Nt))
t_store = zeros(Nt)


i = 0
while t_i < t_max:

	if i % 100 == 0:
		print "t_i =", t_i," |x| =", sqtr(dot(x_i,x_i))
    	#print "x_i =", x_i
    	#print "v_i =", v_i

	# evaluar velocidad relativa
	v_rel = v_f - v_i                      # velocidad relativa
	norm_v_rel = sqtr(dot(v_rel, v_rel))   # norma velocidad relativa
	# evaluar fuerzas sobre la particula
	F_D = 0.5*Cd*norm_v_rel*v_rel          # Fuerza de Drag
	F_i = W                                # Sumatoria de todas las fuerzas que se ejercen sobre la particula

	#print "F_i =", F_i
	# evaluar aceleracion
	a_i = F_i / m

	#print "a_i =", a_i
	# integrar
	x_i_m1 = x_i + v_i*dt + a_i*(dt**2/2)
	v_i_m1 = v_i + a_i*dt

	# avanzar al siguiente paso
	x_store[:, i] = x_i
	v_store[:, i] = v_i
	t_store[:] = t_i

	t_i += dt
	x_i = x_i_m1
	v_i = v_i_m1
	i += 1

# guardar ultimo paso
x_store[:, i] = x_i
v_store[:, i] = v_i
t_store[:] = t_i

figure()
plot(x_store[0, :i], x_store[1, :i])
show()


