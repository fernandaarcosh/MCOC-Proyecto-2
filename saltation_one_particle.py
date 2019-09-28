
# Para simular el salto de una particula, utilizamos la Free flight particle equation, ignorando la fuerza
# de Magnus y la Fuerza historica o de Basset

# Dado que se está simulando el movimiento de la particula en un plano se consideran los ejes X e Y como
# principales


# Se crean listas en donde van a agregagrse l
ejex=[] 
ejey=[]


# 2 situaciones a modelar, la velocidad del flujo junto con la velocidad de la particula

#------------------------------------------
#Formula perfil de velocidad

# Uf(x)= 1/0.41 * ln(30*x)

# A la funcion uno le entrega la ubicacion en el eje x

#--------------------------------------------

#Datos

volumen_particula= 8.80*(10**-9) # [m3]
diametro_particula= 2.1 # [mm]
masa_particula= 2.33*(10**-2) # [gr]

radio_particula= diametro_particula/2

phro= masa_particula/volumen_particula	

coeficiente_drag= 0.47 # se asume particula esferica

coeficiente_lift= (2*radio_particula)/(phro*(volumen_particula**2)*A)

g= 9.8 #[m/s2] ---> gravedad




