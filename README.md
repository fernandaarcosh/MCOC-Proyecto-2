# MCOC-Proyecto-2
Proyecto 2 Métodos Computacionales para Obras Civiles

Integrantes:

- Fernanda Arcos Hernández
- Fabian Cortes Figueroa (https://github.com/fabianszne)
- Roberto Cruz Fernández (https://github.com/RobertoCruzF)
- Anibal Tapia Triviño   (https://github.com/tapiolaa)

# Introducción
En este proyecto se implementará y validará un modelo de simulación de transporte de sedimentos en base a un método lagrangiano, en específico el método de Euler. 
La validación se hará a nivel del comportamiento de una partícula individual y luego con el comportamiento estadístico de cantidades
crecientes de partículas.

# Objetivos
Implementar un modelo de simulacion numerico para transporte de sedimentos de fondo. Comprender aspectos de desempeño de aplicaciones de computación científica tales como IO y complejidad algorítmica.

- [Meta 3] : Implementación y validación del código para una partícula considerando un perfil de velocidades sencillo.
- [Entrega 3] : (script saltation_one_particle.py) Como resultado del código se obtienen dos gráficos en los cuales se observa el movimiento de la partícula, tanto en las direcciones x e y.

- [Meta 4] : Implementacion del código para múltiples partículas considerando un perfil de velocidad sencillo; además de uno complejo.
- [Entrega 4] : (script saltation_many_particles_validation.py) Para el código se comienza por abordar el movimiento de más de una partícula, específicamente 2 con un perfil de velocidad sencillo; luego se pretende poder entregar los resultados para un número de partículas deseadas por el usuario, creando un input y con esto un ciclo que permita guardar las posiciones de cada una de estas partículas.

Fernanda Arcos Hernández:
En la [Entrega 4] como se menciona anteriormente se observan los resultados para el salto de varias particulas en un flujo determinado.A estas partículas se les considera la fuerza de Drag, Lift, Peso y Empuje. Las partículas poseen un diametro de 15 mm; lo que conlleva una densidad de 2650 Kg/m^3 y una velocidad (ustar) predeterminada con los papers entregados de 0.18. 
Se realizó una prueba de como es que funciona el código a medida que aumenta la cantidad de partículas, se observa que el tiempo es mayor a medida que aumenta el número; pero como el código funciona de manera aleatoria, mientras más choques entre particulas existan, mayor sera el tiempo dentro de una misma cantidad; para poder sacar un promedio del tiempo se implemento el código en tres oportunidades por cada n° de particulas, obteniendo un promedio.

En el código se realiza la variación de 1 a 20 párticulas y observar su movimiento durante 5s con un dt =0,001 s.

