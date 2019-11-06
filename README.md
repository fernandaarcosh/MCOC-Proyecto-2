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

En el código se realiza la variación de 1 a 20 párticulas y observar su movimiento durante 2s con un dt =0,001 s.

Gráfico para:

Nparticulas = 2

![Np 2 2s](https://user-images.githubusercontent.com/53495949/66692300-7895fb00-ec73-11e9-88de-44fb537f63e0.JPG)

Nparticulas = 5

![Np5 2s](https://user-images.githubusercontent.com/53495949/66691821-6b770d00-ec6f-11e9-8cc5-5a13a24c4a92.JPG)

Nparticulas = 10 

![Np10 2s](https://user-images.githubusercontent.com/53495949/66691826-7af65600-ec6f-11e9-8653-ed793aa2ff7c.JPG)

Nparticulas = 15

![Np 15 2s](https://user-images.githubusercontent.com/53495949/66691839-8c3f6280-ec6f-11e9-9624-425acd97aa6f.JPG)

Nparticulas = 20


![Np20 2s](https://user-images.githubusercontent.com/53495949/66692314-8c416180-ec73-11e9-9b45-69c4ed135a57.JPG)


El programa para 20 particulas demoro en promedio un total de 414 s

![tiempo](https://user-images.githubusercontent.com/53495949/66692324-b2ff9800-ec73-11e9-8970-d703b92f7979.JPG)

Finalmente las caracteristicas del computador en el cual se ejecuto el código son:
![propiedades equipo](https://user-images.githubusercontent.com/53495949/66691856-bc870100-ec6f-11e9-98ea-c0e38c436041.JPG)

En la [Entrega 6] En esta entrega se pretendió optimizar el código de la [Entrega 4]  
Nuevamente se observa un aumento en el tiempo a medida que aumentan las particulas, pero ahora se ve un comportamiento lineal y una gran mejora en cuanto a la demora del programa; en esta ocasión se hizo correr el programa con 4, 8, 12, 24, 32, 44, 60 y 100 particulas, con un dt de 0.00001 s en un intervalo de 0.05 s . Lo que dio los siguientes resultados.

Np = 32 particulas

![Grafico 32 particulas](https://user-images.githubusercontent.com/53495949/68261447-edfcad80-001e-11ea-8e51-ffe84fe4b1bd.png)

Np = 60 particulas

![Grafico 60 particulas](https://user-images.githubusercontent.com/53495949/68261582-56e42580-001f-11ea-96df-2f1a725d23ac.png)

Np = 100 particulas

![Grafico 100 particulas](https://user-images.githubusercontent.com/53495949/68261922-631cb280-0020-11ea-9823-80a7da754e43.png)


Finalmente se obtiene la siguiente relación:

Gráfico de tiempo
![Captura](https://user-images.githubusercontent.com/53495949/68262163-2f8e5800-0021-11ea-8a31-305994a38e42.JPG)

