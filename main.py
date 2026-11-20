# Solucion de la ecuacion de Blasius mediante RK4

import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from numpy import array, arange, floor, ceil
import pandas as pd
import os

# Parametros de la solucion para la ecuacion de Blasius

h_list = [0.1, 0.05, 0.01, 0.001] # Tamano de paso
eta_max = [6, 8, 10, 12] # Ultimo valor de eta 

max_iter = 40 # Numero maximo de iteraciones para el metodo Newton-Raphson
delta = 1e-8 # Incremento para la derivada numerica 
tol = 1e-12 # Tolerancia para Newton-Raphson 

# Parametros de la solucion para la ecuacion termica

h_temp = 0.01   # Tamano de paso (rapido y preciso)
max_iter_temp = 20

# Metodo Newton-Raphson y guardado para t'(0)

h_pr = 0.01 # Tamano de paso para el numero de Prandtl
Pr_max = 10  # Numero de Prandtl maximo
Pr_max_perfil = 1 # Numero de Prandtl maximo para graficar perfiles
Pr = arange(h_pr, Pr_max + h_pr, h_pr) # Empieza en 0.01 para evitar Pr = 0
tp = 0.2 # Semilla para t'(0)
tp_list = []

# Direccion de guardado de datos
ruta = os.path.join(os.getcwd(), "resultados") # type: ignore
os.makedirs(ruta, exist_ok=True)

# Integracion Blasius

def RK4(fpp0, h, eta_max):
    eta = 0
    f, fp, fpp = 0.0, 0.0, fpp0
    
    eta_var = [eta]
    f_var = [f]
    fp_var = [fp]
    fpp_var = [fpp]
    
    N = int(round(eta_max / h))
    h2 = 0.5 * h
    h6 = h / 6.0
    
    for _ in range(N):
        # k1
        k1_0 = fp
        k1_1 = fpp
        k1_2 = -0.5 * f * fpp
        
        # k2
        f2 = f + h2 * k1_0
        fp2 = fp + h2 * k1_1
        fpp2 = fpp + h2 * k1_2
        k2_0 = fp2
        k2_1 = fpp2
        k2_2 = -0.5 * f2 * fpp2
        
        # k3
        f3 = f + h2 * k2_0
        fp3 = fp + h2 * k2_1
        fpp3 = fpp + h2 * k2_2
        k3_0 = fp3
        k3_1 = fpp3
        k3_2 = -0.5 * f3 * fpp3
        
        # k4
        f4 = f + h * k3_0
        fp4 = fp + h * k3_1
        fpp4 = fpp + h * k3_2
        k4_0 = fp4
        k4_1 = fpp4
        k4_2 = -0.5 * f4 * fpp4
        
        f += h6 * (k1_0 + 2.0 * k2_0 + 2.0 * k3_0 + k4_0)
        fp += h6 * (k1_1 + 2.0 * k2_1 + 2.0 * k3_1 + k4_1)
        fpp += h6 * (k1_2 + 2.0 * k2_2 + 2.0 * k3_2 + k4_2)
        eta += h
        
        eta_var.append(eta)
        f_var.append(f)
        fp_var.append(fp)
        fpp_var.append(fpp)
        
    return [eta_var, f_var, fp_var, fpp_var]

# Encontrar la f''(0) que cumpla con f'(infinito) = 1

def R(fpp0, h, eta_max):
    _, _, fp, _ = RK4(fpp0, h, eta_max)
    return fp[-1] - 1.0

# Metodo Newton-Rapshon para distintas combinaciones de h y eta_max

min_error = 1 # Sirve para encontrar la condicion f''(0) de minimo error
opt_fpp = []  # Lista para almacenar la mejor combinacion de h y eta_max, basados en la f''(0) con el minimo error 

with open(os.path.join(ruta, 'fpp.txt'), 'w', encoding='UTF-8') as save: # type: ignore
    save.write("Resultados finales para la tercer condicion inicial\n\n")

    for i in h_list:
        
        for j in eta_max:
            
            fpp = 0.1 # Semilla para f''(0)

            for k in range(max_iter):
                gpp = R(fpp, i, j)
                
                if abs(gpp) <= tol:
                    save.write(f"h = {i:<6g} | eta_max = {j:<3d} | iteraciones NR = {k:<2d} | ")
                    break
                
                # Derivada numerica (Diferencia central)
                dg = (R(fpp + delta, i, j) - R(fpp - delta, i, j))/(2*delta)

                # Newton  
                fpp = fpp - gpp / dg
                error = abs(gpp)
                
            # Guarda la combinacion de menor error
            if error < min_error:
                min_error = error
                opt_fpp = [float(fpp), i, j]
                
            save.write(f"f''(0) = {fpp:<19.16f} | Error = {error}\n")


# Guardar perfil completo de la mejor solucion
eta_var, f_var, fp_var, fpp_var = RK4(opt_fpp[0], opt_fpp[1], opt_fpp[2])

with open(os.path.join(ruta, 'Blasius_solution.txt'), 'w', encoding='UTF-8') as save: # type: ignore
    save.write(f"{'eta':<8} | {'f':<12} | {'f\'':<12} | {'f\'\'':<12}\n")
    save.write("-" * 52 + "\n")
    for eta_i, f_i, fp_i, fpp_i in zip(eta_var, f_var, fp_var, fpp_var):
        save.write(f"{eta_i:<8.3f} | {f_i:<12.6f} | {fp_i:<12.6f} | {fpp_i:<12.6f}\n")


# Valor de eta maximo para los perfiles de temperatura

eta_max_temp = opt_fpp[2]  # Numero de pasos (eta_max = 10)

# Generacion de graficas

plt.rcParams.update({
    'font.size': 11,
    'axes.labelsize': 12,
    'axes.titlesize': 13,
    'legend.fontsize': 10,
    'xtick.direction': 'in',
    'ytick.direction': 'in'
})

# 1. Grafica combinada (f, f' y f'' vs eta)

fig, ax = plt.subplots(figsize=(8, 5.5), dpi=300)

ax.plot(eta_var, f_var, label=r"$f(\eta)$ (Función de corriente)", color='#1f77b4', lw=2)
ax.plot(eta_var, fp_var, label=r"$f'(\eta) = u/U_\infty$ (Velocidad)", color='#d62728', lw=2)
ax.plot(eta_var, fpp_var, label=r"$f''(\eta)$ (Esfuerzo cortante)", color='#2ca02c', lw=2)

ax.axhline(1.0, color='gray', linestyle='--', linewidth=1, alpha=0.7)

ax.set_title(rf"Solución de Blasius ($h = {opt_fpp[1]}$, $\eta_{{\max}} = {opt_fpp[2]}$, $f''(0) = {opt_fpp[0]:.5f}$)")
ax.set_xlabel(r"Variable de similitud, $\eta$")
ax.set_ylabel(r"$f(\eta),\ f'(\eta),\ f''(\eta)$")
ax.set_xlim(0, opt_fpp[2])
ax.set_ylim(0, max(f_var) * 1.05)

ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.grid(True, which='major', linestyle='-', alpha=0.4)
ax.grid(True, which='minor', linestyle=':', alpha=0.2)
ax.legend(frameon=True, facecolor='white', framealpha=0.9)

plt.tight_layout()
plt.savefig(os.path.join(ruta, 'Blasius_combinada.png'), dpi=300) # type: ignore
plt.close(fig)

# 2. Tres graficas individuales (f vs eta, f' vs eta, f'' vs eta)

graficas = [
    (f_var,   r"$f(\eta)$",   r"Función de corriente adimensional, $f(\eta)$", '#1f77b4', 'Blasius_f_vs_eta.png'),
    (fp_var,  r"$f'(\eta)$",  r"Perfil de velocidad adimensional, $f'(\eta) = u/U_\infty$", '#d62728', 'Blasius_fp_vs_eta.png'),
    (fpp_var, r"$f''(\eta)$", r"Derivada segunda (esfuerzo cortante), $f''(\eta)$", '#2ca02c', 'Blasius_fpp_vs_eta.png')
]

for datos_y, ylabel, titulo, color, nombre_archivo in graficas:
    fig, ax = plt.subplots(figsize=(7, 5), dpi=300)
    
    ax.plot(eta_var, datos_y, color=color, lw=2.2)
    
    # Marcar el espesor de capa limite (eta donde f' = 0.99) en la grafica de f'
    if datos_y is fp_var:
        ax.axhline(0.99, color='black', linestyle='--', lw=1, alpha=0.7, label=r"$u/U_\infty = 0.99$ ($\eta \approx 4.91$)")
        ax.axvline(4.91, color='black', linestyle=':', lw=1, alpha=0.7)
        ax.legend(loc='lower right')

    ax.set_title(titulo)
    ax.set_xlabel(r"Variable de similitud, $\eta$")
    ax.set_ylabel(ylabel)
    ax.set_xlim(0, opt_fpp[2])
    ax.set_ylim(0, max(datos_y) * 1.08)
    
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.grid(True, which='major', linestyle='-', alpha=0.4)
    ax.grid(True, which='minor', linestyle=':', alpha=0.2)
    
    plt.tight_layout()
    plt.savefig(os.path.join(ruta, nombre_archivo), dpi=300) # type: ignore
    plt.close(fig)
    
# Solucion a ecuacion de energia

def G(y, f, Pr):
    return [y[1], -0.5*Pr*y[1]*f]

def paso_Temp(y, h, f, Pr):
    y = array(y)
    k1 = array(G(y, f, Pr))
    k2 = array(G(y + (h/2)*k1, f, Pr))
    k3 = array(G(y + (h/2)*k2, f, Pr))
    k4 = array(G(y + h*k3, f, Pr))
    return y + (h/6)*(k1 + 2*k2 + 2*k3 + k4)

eta_var, f_var, fp_var, fpp_var = RK4(opt_fpp[0], h_temp, eta_max_temp)

def RK4_Temp(y, h, f, Pr, eta_max):
    t_var, tp_var = [y[0]], [y[1]]
    N = int(round(eta_max/h))
    for i in range(N):
        y = paso_Temp(y, h, f[i], Pr)
        t_var.append(y[0])
        tp_var.append(y[1])
    return [t_var, tp_var]

def g(tp0, f, Pr):
    y0 = [0, tp0]
    t, _ = RK4_Temp(y0, h_temp, f, Pr, eta_max_temp)
    return t[-1] - 1.0


with open(os.path.join(ruta, 'tp.txt'), 'w', encoding='UTF-8') as save: # type: ignore
    save.write(f"{'Pr':<6} | {'Iteraciones NR':<14} | {'t\'(0)':<20} | {'Error':<14}\n")
    save.write("-" * 64 + "\n")

    for pr in Pr:
        for n in range(max_iter_temp):
            gp = g(tp, f_var, pr)
            error = abs(gp)

            if error < tol:
                tp_list.append(tp)
                break

            dt = (g(tp + delta, f_var, pr) - g(tp - delta, f_var, pr)) / (2*delta)
            tp = tp - gp / dt

        save.write(f"{pr:<6.2f} | {n:<14d} | {tp:<20.16f} | {error:<14.6e}\n")

# Grafica comparativa de Pohlhausen
dtp = 0.332 * (Pr)**(1/3)

plt.plot(Pr, tp_list, color='green', label="Resultado numérico")
plt.plot(Pr, dtp, color='red', linestyle='--', label=r"Aproximación de Pohlhausen ($0.332\,Pr^{1/3}$)")
plt.title(r"Aproximación vs Resultado numérico de $\theta'(0)$")
plt.xlabel("Número de Prandtl, Pr")
plt.ylabel(r"$\theta'(0)$")
plt.xlim(0, Pr_max)
plt.ylim(0)
plt.grid(True, color="gray", alpha=0.4)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(ruta, 'Pr_vs_tp0.png'), dpi=300) # type: ignore
plt.close()

fig, ax = plt.subplots(figsize=(8, 5.5), dpi=300)

for idx_tp, tp_val in enumerate(tp_list):
    # Graficar cada 10 pasos: Pr = 0.10, 0.20, ..., Pr_max (y el primero Pr = 0.01)
    if idx_tp == 0 or (idx_tp + 1) % 10 == 0 and Pr[idx_tp] <= Pr_max_perfil:
        t_var, tp_var = RK4_Temp([0, tp_val], h_temp, f_var, Pr[idx_tp], eta_max_temp)
        ax.plot(eta_var, t_var, lw=1.8, label=rf"$Pr = {Pr[idx_tp]:.2f}$")

ax.axhline(1.0, color='gray', linestyle='--', lw=1, alpha=0.7)
ax.set_title(r"Perfiles de temperatura adimensional $\theta(\eta)$ para distintos números de Prandtl")
ax.set_xlabel(r"Variable de similitud, $\eta$")
ax.set_ylabel(r"Temperatura adimensional, $\theta(\eta)$")
ax.set_xlim(0, eta_max_temp)
ax.set_ylim(0, 1.05)
ax.grid(True, linestyle='-', alpha=0.4)
ax.legend(loc='lower right', ncol=2, fontsize=9)

plt.tight_layout()
plt.savefig(os.path.join(ruta, 'Perfiles_Temperatura_Pr.png'), dpi=300) # type: ignore
plt.close(fig)

# Graficas de perfil de capa limite (Estilo Fisico) para Pr = {0.6, 1.0, 5.0}

for pr_obj in [0.6, 1.0, 5.0]:
    # Buscar theta'(0) para este Pr
    tp_obj = 0.332 * (pr_obj ** (1/3))
    for _ in range(max_iter_temp):
        gp = g(tp_obj, f_var, pr_obj)
        if abs(gp) < tol:
            break
        dt = (g(tp_obj + delta, f_var, pr_obj) - g(tp_obj - delta, f_var, pr_obj)) / (2 * delta)
        tp_obj = tp_obj - gp / dt

    t_var, _ = RK4_Temp([0, tp_obj], h_temp, f_var, pr_obj, eta_max_temp)

    fig, ax = plt.subplots(figsize=(8, 6), dpi=300)

    ax.plot(eta_var, fp_var, color='black', linestyle='-', lw=1.8,
            label=r"Velocidad $u/u_\infty$")
    ax.plot(eta_var, t_var, color='red', linestyle='--', lw=1.8,
            label=rf"Temperatura $\theta$ (Pr={pr_obj})")

    ax.text(2.20, 0.52, f"f''(0) = {opt_fpp[0]:.5f}", color='black', fontsize=10)
    ax.text(2.20, 0.44, rf"$\theta'(0)$ = {tp_obj:.5f}", color='red', fontsize=10)

    ax.set_title("Perfil de Capa Límite", fontsize=12)
    ax.set_xlabel(r"Variable de similitud $\eta$ (Distancia desde la pared)", fontsize=11)
    ax.set_ylabel("Magnitud Adimensional", fontsize=11)
    ax.set_xlim(0, 8)
    ax.set_ylim(0.0, 1.1)
    ax.grid(True, linestyle=':', alpha=0.5)
    ax.legend(loc='lower right', frameon=True, fontsize=10)

    plt.tight_layout()
    nombre_archivo = f"Perfil_Capa_Limite_Pr_{str(pr_obj).replace('.', '_')}.png"
    plt.savefig(os.path.join(ruta, nombre_archivo), dpi=300) # type: ignore
    plt.close(fig)

print("Graficas guardadas exitosamente en:", ruta)