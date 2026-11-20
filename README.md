# blasius-boundary-layer-rk4
Numerical similarity solution for velocity and thermal boundary layers using the shooting method and RK4 (RMF-E).
*[Read this documentation in English](README_EN.md)*

# Solución Numérica de Similitud para las Capas Límite de Velocidad y Térmica sobre una Placa Plana

Este repositorio contiene el código fuente unificado, los datos de salida tabulados y las rutinas de visualización asociados al manuscrito enviado al *Suplemento de la Revista Mexicana de Física (RMF-E)*.

## Descripción del Modelo

El programa resuelve el problema clásico acoplado de la capa límite laminar, bidimensional, incompresible y en régimen permanente sobre una placa plana isotérmica mediante transformaciones de similitud (ecuación de Blasius y ecuación de energía térmica):

1. **Capa límite hidrodinámica:** Reducción de la ecuación de Blasius ($f''' + \frac{1}{2}ff'' = 0$) a un sistema autónomo de primer orden, integrado mediante un esquema de **Runge-Kutta de cuarto orden (RK4)** de paso uniforme y acoplado al método de disparo (*Shooting Method*) con ajuste de raíces mediante **Newton-Raphson** (derivadas aproximadas por diferencias centrales).
2. **Análisis de sensibilidad espacial y asintótica:** Evaluación sistemática de 16 combinaciones de tamaño de paso $h \in \{0.1, 0.05, 0.01, 0.001\}$ y dominio asintótico $\eta_{\max} \in \{6, 8, 10, 12\}$ para identificar la configuración óptima ($h = 0.01$, $\eta_{\max} = 8$, $|R(s)| = 1.15 \times 10^{-12}$).
3. **Capa límite térmica:** Resolución acoplada de la ecuación de energía ($\theta'' + \frac{Pr}{2}f\theta' = 0$) para el espectro de números de Prandtl $0.01 \le Pr \le 10.0$, validando los gradientes térmicos superficiales $\theta'(0)$ frente a valores de referencia y la correlación analítica de Pohlhausen ($0.332\,Pr^{1/3}$).

---

## Requisitos e Instalación

* **Versión de Python recomendada:** Python 3.10 o superior (probado en Python 3.11).
* **Dependencias:** `numpy`, `matplotlib`.

### Instalación rápida

Clona este repositorio e instala las dependencias necesarias ejecutando:

```bash
git clone https://github.com/Scienceeeee/blasius-boundary-layer-rk4.git
cd blasius-boundary-layer-rk4
pip install -r requirements.txt
