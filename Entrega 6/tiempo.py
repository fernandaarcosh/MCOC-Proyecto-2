from matplotlib.pylab import *
from matplotlib import pyplot

x = [4, 8, 12, 24, 32, 44, 60, 100] #Nparticulas
y = [11.817, 18.632, 27.036, 56.889, 85.397, 114.753, 169.319, 314.969 ] # Tiempo que se demora total

pyplot.plot(x,y, "--")
pyplot.plot(x,y, "o")
plt.xlabel("Numero Particulas") # Unidades eje x
plt.ylabel("Tiempo [s]")	# Unidades eje y

plt.title("Tiempo vs N particulas")
show()