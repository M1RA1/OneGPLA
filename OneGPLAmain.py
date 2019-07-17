import numpy as num
import OneGPLAfunctions as PIC

#DINAMICA == 0 no se grafica al correr, DINAMICA == 1 grafica
DINAMICA=0
#Esta función se utiliza para transportar con facilidad 
#el conjunto de variables 
def part(x, v, tag, q, m):
    return (x, v, tag, q, m, q / m)
#Cargar las semillas creadas en OneGPLA semilla
x1 = num.load('Semilla/semillaposicion1.npy')
v1 = num.load('Semilla/semillavelocidad1.npy')
x2 = num.load('Semilla/semillaposicion2.npy')
v2 = num.load('Semilla/semillavelocidad2.npy')
xi = num.load('Semilla/semillaposicioni.npy')
vi = num.load('Semilla/semillavelocidadi.npy')
numtime, moveions, Npb1,Npb2,Npi,NpT,L,vd1,vd2,vt1,vt2,vti,ni,q1,q2,qi,m1,m2,mi,EPS0,dh,dt,vmag0, vmag1, vmag2= numpy.load('Semilla/semilladatos.npy')
numtime = int(numtime)
Npb1 = int(Npb1)
Npb2 = int(Npb2)
Npi = int(Npi)
q0 = qi
ni = int(ni)
np = int(Npb1)
parts = []  # vector vacío para los electrones
partsi = []  # vector vacío para los iones

for p in range(0, Npb1):
    parts.append(part(x1[p], v1[p], 0, q1, m1))
for p in range(0, Npb2):
    parts.append(part(x2[p], v2[p], 1, q2, m2))
for p in range(0, Npi):
    partsi.append(part(xi[p], vi[p], 0, qi, mi))

dh = L / (ni - 1)  # Tamaño de la celda

# Resetear los campos
l2 = 0
phi = []
ef = []
rho = []
rhoi = []
rho2 = []
x = []

for i in range(0, ni):
    phi.append(0)
    ef.append(0)
    rho.append(0)
    rhoi.append(0)
    rho2.append(0)
    x.append(i * dh)
 
it = 0         # Contador
E_acumulador = numpy.empty((ni, numtime + 1))
phi_acumulador = numpy.empty((ni, numtime + 1))
rho_acumulador = numpy.empty((ni, numtime + 1))
EnergiaK = []
EnergiaP = []
EnergiaT = []
VT = []
for it in range(0, numtime + 1):
    if it % 10 == 0:
        print('it', it)
    if moveions==0:
        rho, phi, ef, parts = PIC.stepion(ni, rho, rhoi, phi, ef, dt, dh, parts, partsi, EPS0, L, it)
    if moveions==1:
        rho, phi, ef, parts, partsi = PIC.stepmoveion(ni, rho, rhoi, phi, ef, dt, dh, parts, partsi, EPS0, L, it)
    for uni in range(0, ni):
        E_acumulador[uni, it] = ef[uni]
        phi_acumulador[uni,it] = phi[uni]
        rho_acumulador[uni, it] = rho[uni]
    K, P = PIC.energias(rho, phi, parts)
    EnergiaK.append(K)
    EnergiaP.append(P)
    EnergiaT.append(K+P)
    VT.append(it)
    #Gráfica Dinámica
    if it % 5 == 0 and DINAMICA == 1:
        PIC.plots(parts, L, vmag0)        
    #Se guardan las posiciones cada 25 pasos
    if it % 25 == 0:
        xa1 = []
        xa2 = []
        xai = []
        va1 = []
        va2 = []
        vai = []
        for p in range(0, (Npb1 + Npb2)):
            P = parts[p]
            x, v, tag, q, m, qm = P
            if tag ==0:
                xa1.append(x)
                va1.append(v)
            else:
                xa2.append(x)
                va2.append(v)
        for p in range(0, Npi):
            Pi = partsi[p]
            x, v, tagi, qi, mi, qmi = Pi
            xai.append(x)
            vai.append(v)
        num.save('Posiciones/1posicion{0}.npy'.format(it), xa1)
        num.save('Velocidades/1velocidad{0}.npy'.format(it), va1)
        num.save('Posiciones/2posicion{0}.npy'.format(it), xa2)
        num.save('Velocidades/2velocidad{0}.npy'.format(it), va2)
        num.save('Posiciones/posicionion{0}.npy'.format(it), xai)
        num.save('Velocidades/velocidadion{0}.npy'.format(it), vai)
        num.save('Campo/campo{0}.npy'.format(it), E_acumulador)
        num.save('Potencial/phi{0}.npy'.format(it), phi_acumulador)
        num.save('Densidad/rho{0}.npy'.format(it), rho_acumulador)
        num.save('Energías/EnergiaK{0}.npy'.format(it), EnergiaK)
        num.save('Energías/EnergiaP{0}.npy'.format(it), EnergiaP)





