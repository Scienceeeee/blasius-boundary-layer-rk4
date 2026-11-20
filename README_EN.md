# blasius-boundary-layer-rk4
Numerical similarity solution for velocity and thermal boundary layers using the shooting method and RK4 (RMF-E).

*[Leer esta documentación en Español](README.md)*

# Numerical Similarity Solution for Velocity and Thermal Boundary Layers on a Flat Plate

This repository contains the unified source code, tabulated output data, and visualization routines associated with the manuscript submitted to *Suplemento de la Revista Mexicana de Física (RMF-E)*.

## Model Description

The code solves the classical coupled problem of a steady, two-dimensional, incompressible laminar boundary layer over an isothermal flat plate using similarity transformations (Blasius equation and thermal energy equation):

1. **Hydrodynamic boundary layer:** Reduction of the Blasius equation ($f''' + \frac{1}{2}ff'' = 0$) to an autonomous first-order system, integrated via a uniform-step **fourth-order Runge-Kutta (RK4)** scheme coupled with the **Shooting Method** and root-finding via **Newton-Raphson** (derivatives approximated using central finite differences).
2. **Spatial and asymptotic sensitivity analysis:** Systematic evaluation of 16 combinations of step size $h \in \{0.1, 0.05, 0.01, 0.001\}$ and asymptotic domain $\eta_{\max} \in \{6, 8, 10, 12\}$ to determine the optimal discretization ($h = 0.01$, $\eta_{\max} = 8$, $|R(s)| = 1.15 \times 10^{-12}$).
3. **Thermal boundary layer:** Coupled solution of the thermal energy equation ($\theta'' + \frac{Pr}{2}f\theta' = 0$) across the Prandtl number spectrum $0.01 \le Pr \le 10.0$, validating wall temperature gradients $\theta'(0)$ against benchmark literature values and the analytical Pohlhausen correlation ($0.332\,Pr^{1/3}$).

---

## Requirements & Installation

* **Recommended Python version:** Python 3.10 or higher (tested on Python 3.11).
* **Dependencies:** `numpy`, `matplotlib`.

### Quick Installation

Clone this repository and install the required packages by running:

```bash
git clone https://github.com/Scienceeeee/blasius-boundary-layer-rk4.git
cd blasius-boundary-layer-rk4
pip install -r requirements.txt
