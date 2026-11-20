================================================================================

SOLUCIÓN NUMÉRICA DE SIMILITUD PARA LAS CAPAS LÍMITE DE VELOCIDAD Y TÉRMICA

SOBRE UNA PLACA PLANA MEDIANTE EL MÉTODO DE DISPARO Y RUNGE-KUTTA DE CUARTO ORDEN



NUMERICAL SIMILARITY SOLUTION FOR VELOCITY AND THERMAL BOUNDARY LAYERS

ON A FLAT PLATE VIA THE SHOOTING METHOD AND FOURTH-ORDER RUNGE-KUTTA

================================================================================



Autores / Authors:

- Abed Gael Hinojosa Baez

- Germán Reza Bass

- Edoardo García Díaz



--------------------------------------------------------------------------------

[ESPAÑOL] INSTRUCCIONES DE USO

--------------------------------------------------------------------------------



1. DESCRIPCIÓN GENERAL

Este repositorio contiene el código fuente en Python para resolver de forma numérica

las ecuaciones diferenciales ordinarias acopladas de Blasius (capa límite de 

velocidad) y de energía (capa límite térmica) sobre una placa plana isotérmica, 

utilizando el método de disparo (Shooting Method) combinado con el esquema de 

integración de Runge-Kutta de cuarto orden (RK4) y el método de Newton-Raphson.



2. REQUISITOS E INSTALACIÓN

El código requiere Python 3.x y las bibliotecas científicas estándar:

- numpy

- matplotlib



Si no las tienes instaladas, puedes instalarlas ejecutando en tu terminal:

    pip install numpy matplotlib



3. EJECUCIÓN

El script está listo para ejecutarse directamente sin configuración adicional:

    python nombre_del_script.py



Al ejecutarse, el programa calculará automáticamente las soluciones de similitud, 

convergerá las condiciones iniciales faltantes en la pared (f''(0) y t'(0)) y 

generará las gráficas de los perfiles de velocidad y temperatura correspondientes.



4. PARÁMETROS CONFIGURABLES DE LA SIMULACIÓN

Si deseas experimentar con la precisión numérica, el dominio o los rangos físicos, 

puedes modificar las siguientes variables directamente en el bloque de parámetros 

del código:



A) Parámetros para la ecuación de Blasius (Capa límite de velocidad):

- h_list = [0.1, 0.05, 0.01, 0.001]

  Lista de tamaños de paso (Delta eta) utilizados para evaluar la convergencia 

  y precisión del método RK4.

- eta_max = [6, 8, 10, 12]

  Lista de valores máximos del dominio de similitud (eta -> infinito) donde se 

  impone la condición de frontera de corriente libre (f'(eta_max) = 1).

- max_iter = 40

  Número máximo de iteraciones permitidas para el algoritmo de Newton-Raphson.

- delta = 1e-8

  Incremento utilizado para aproximar la derivada numérica mediante diferencias 

  finitas dentro del método de disparo.

- tol = 1e-12

  Tolerancia de error absoluto requerida para considerar que Newton-Raphson ha 

  convergido.



B) Parámetros para la ecuación térmica (Capa límite térmica):

- h_temp = 0.01

  Tamaño de paso (Delta eta) empleado para integrar la ecuación de energía 

  (compromiso óptimo entre rapidez y precisión).

- max_iter_temp = 20

  Número máximo de iteraciones de Newton-Raphson para ajustar la pendiente térmica 

  en la pared t'(0).



C) Parámetros del número de Prandtl (Pr) y graficación:

- h_pr = 0.01

  Incremento o tamaño de paso para barrer el rango de números de Prandtl.

- Pr_max = 10

  Valor máximo del número de Prandtl hasta el cual se calcula y guarda t'(0).

- Pr_max_perfil = 1

  Valor máximo del número de Prandtl utilizado para graficar los perfiles 

  espaciales de temperatura.

- Pr = arange(h_pr, Pr_max + h_pr, h_pr)

  Arreglo de valores de Prandtl evaluados (inicia en h_pr = 0.01 para evitar la 

  singularidad en Pr = 0).

- tp = 0.2

  Valor semilla (estimación inicial) para la derivada de la temperatura 

  adimensional en la pared, t'(0), al iniciar el método de disparo.





--------------------------------------------------------------------------------

[ENGLISH] USAGE INSTRUCTIONS

--------------------------------------------------------------------------------



1. OVERVIEW

This repository contains the Python source code to numerically solve the coupled 

ordinary differential equations for Blasius (velocity boundary layer) and energy 

(thermal boundary layer) over an isothermal flat plate, using the Shooting Method 

combined with the fourth-order Runge-Kutta (RK4) integration scheme and the 

Newton-Raphson root-finding algorithm.



2. REQUIREMENTS AND INSTALLATION

The script requires Python 3.x and standard scientific libraries:

- numpy

- matplotlib



If not installed, you can install them via terminal:

    pip install numpy matplotlib



3. EXECUTION

The script is ready to run out-of-the-box without further setup:

    python script_name.py



Upon execution, the program automatically computes the similarity solutions, 

converges the missing wall boundary conditions (f''(0) and t'(0)), and generates 

the plots for the velocity and temperature profiles.



4. CONFIGURABLE SIMULATION PARAMETERS

If you wish to explore numerical accuracy, domain size, or physical regimes, you 

can modify the following variables directly in the parameter section of the code:



A) Blasius Equation Parameters (Velocity boundary layer):

- h_list = [0.1, 0.05, 0.01, 0.001]

  List of step sizes (Delta eta) used to test the convergence and accuracy of 

  the RK4 integrator.

- eta_max = [6, 8, 10, 12]

  List of upper bounds for the similarity domain (eta -> infinity) where the 

  free-stream boundary condition (f'(eta_max) = 1) is enforced.

- max_iter = 40

  Maximum number of iterations allowed for the Newton-Raphson solver.

- delta = 1e-8

  Perturbation increment used to approximate the numerical derivative via finite 

  differences in the shooting method.

- tol = 1e-12

  Absolute error tolerance required for Newton-Raphson convergence.



B) Thermal Equation Parameters (Thermal boundary layer):

- h_temp = 0.01

  Step size (Delta eta) used to integrate the energy equation (selected for a 

  balance between speed and accuracy).

- max_iter_temp = 20

  Maximum number of Newton-Raphson iterations to find the wall temperature 

  gradient t'(0).



C) Prandtl Number (Pr) and Plotting Parameters:

- h_pr = 0.01

  Step increment used to sweep across the Prandtl number range.

- Pr_max = 10

  Maximum Prandtl number value for computing and storing t'(0).

- Pr_max_perfil = 1

  Maximum Prandtl number value included when plotting temperature profiles.

- Pr = arange(h_pr, Pr_max + h_pr, h_pr)

  Array of evaluated Prandtl numbers (starts at h_pr = 0.01 to avoid the 

  singularity at Pr = 0).

- tp = 0.2

  Initial guess (seed) for the dimensionless wall temperature gradient, t'(0), 

  to initialize the shooting method.

================================================================================
