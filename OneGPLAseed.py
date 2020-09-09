import matplotlib.pyplot as plt
import numpy as np
import OneGPLAfunctions as PIC

# Número de partículas de los haces
numtime = 300 # Número de iteraciones a ejecutar
moveions = 0  # 0 para que no se muevan, 1 para que se muevan.
Npb1 = 2000 # Partíulas beam 1
Npb2 = 2000 # Partíulas beam 2
Npi = 2000  # Iones
NpT = Npb1+Npb2 + Npi #Total partíclas
L = 2* np.pi #Longitud rejilla computacional
# Densidad lineal de partículas de los haces
n1 = Npb1/ L 
n2 = Npb2/ L
ni = Npi / L
# Velocidades de desfase de los haces
vd1= -1
vd2= 1
vt1= 0.1 # El 1 es para las velocidades paralelas y el 2 para las perpendiculares
vt2= 0.1
vti= 0.01 # Velocidad térmica iónica

num_cells = 64 # Número de celdas
Nicell = num_cells + 1
dt = 0.05  #Paso temporal
dh = L / (ni - 1)  # Espacio entre celdas
qm1 = -1   # Carga entre masa de las partículas 
qm2 = -1  
qmi = 1/100
# Frecuencia del plasma
wp = 1  
EPS0 = 1 # Se normaliza la permitividad eléctrica a 1 para evitar overflow
m1=1 #Masa particulas beam 1
m2=1 #Masa particulas beam 2
mi = 100 #Masa iones

# Carga calculada a partir de la frecuencia del plasma
q1 =  wp**2 * (1 / qm1) * EPS0 * L / (NpbT)
q2 = q1
# charge of the neutralizing background
qi = -0.5 * (q1 + q2)

# masses, only used for energy

vmag1 = 1  # de v1 para graficar
vmag2 = -1  # de V2 para graficar
vmag0 = 0.5 * (np.abs(vmag1) + np.abs(vmag2)) # vmag0 = 1???

#Crea estructura con todos los parámetros para almacenarlos.
semilladatos= numtime, moveions, Npb1, Npb2, Npi, NpT, L, vd1, vd2, vt1, vt2, vti, Nicell, q1, q2, qi, m1, m2, mi, EPS0, dh, dt,vmag0, vmag1, vmag2
#En esta parte se crea la semilla
positions1 = PIC.CargarRandomPosicion(Npb1, L)
positions2 = PIC.CargarRandomPosicion(Npb2, L)
positionsi = PIC.CargarRandomPosicion(Npi, L)
velocities1 = PIC.CargarVelocidadMaxwell1D(Npb1, n1, vd1, vt1)
velocities2 = PIC.CargarVelocidadMaxwell1D(Npb2, n2, vd2, vt2)
velocitiesi = PIC.CargarVelocidadMaxwell1D(Npi, ni, 0, vti)
positions1 = PIC.CondicionesDeFrontera(Npb1, positions1, L)
positions2 = PIC.CondicionesDeFrontera(Npb2, positions2, L)
positionsi = PIC.CondicionesDeFrontera(Npi, positionsi, L)
positions0 = np.append(positions1,positions2)
velocities0 = np.append(velocities1,velocities2)

#Repartimos uniformemente los iones
delta_x = L / Npi
mode = 1  # de formulario
dx0 = 0.0001  # de formulario
for p in range(0, Npi):
    x0 = (p + 0.5) * delta_x
    # perturbed positions
    theta = 2 * np.pi * mode * x0 / L
    dx = dx0 * np.cos(theta)
    x1 = x0 + dx
    if x1 < 0:
        x1 = x1 + L
    if x1 >= L:
        x1 = x1 - L
        
    positionsi[p] = x1

#!
plt.scatter(positions1, velocities1)
plt.scatter(positions2, velocities2)
plt.scatter(positionsi, velocitiesi)

#Se guardan las semillas en la carpeta Semilla
np.save('Semilla/semillaposicion.npy', positions0)
np.save('Semilla/semillavelocidad1.npy', velocities0)
np.save('Semilla/semillaposicion1.npy', positions1)
np.save('Semilla/semillavelocidad1.npy', velocities1)
np.save('Semilla/semillaposicion2.npy', positions2)
np.save('Semilla/semillavelocidad2.npy', velocities2)
np.save('Semilla/semillaposicioni.npy', positionsi)
np.save('Semilla/semillavelocidadi.npy', velocitiesi)
np.save('Semilla/semilladatos.npy', semilladatos)
