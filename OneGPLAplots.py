import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate
from mpl_toolkits.mplot3d import Axes3D
import math
import Dispersion as D
from scipy.fftpack import fft2              # FFT dos dimensional
from scipy.fftpack import fftshift as shift # Corrimiento al cero 
import OneGPLAfunctions as OD


tfield= 100 # graficar las magnitudes en intervalos de tfield 
tiempo=300   #Último archivo a graficar
dt=0.025
L = 32*math.pi
N = 64
dh = L/N
npe=533688
Np = 743403
Npvel = 100  # Numero de puntos de velocidad X
v0 = 10
dv = 2 * v0 / Npvel
ve = sp.arange(-v0, v0, dv)
X=np.arange(0,L+dh,dh)
T=dt * np.arange(0,tiempo+1)

phi = np.load('Potencial/phi{0}.npy'.format(tiempo))
campo = np.load('Campo/campo{0}.npy'.format(tiempo))
rho = np.load('Densidad/rho{0}.npy'.format(tiempo))
EK = np.load('Energías/EnergiaK{0}.npy'.format(tiempo))
EP = np.load('Energías/EnergiaP{0}.npy'.format(tiempo))
#Se guardan las graficas de los espacios de fase 
# y funciones de distribución a razón de it % 100 

#GRAFICA DEL ESPACIO DE FASE Y LA FUNCIÓN DE DISTRIBUCIÓN
for it in range(0,tiempo+1):
    
    if it % 100 == 0:
       x1 = np.load('Posiciones/1posicion{0}.npy'.format(it))
       v1 = np.load('Velocidades/1velocidad{0}.npy'.format(it))
       x2 = np.load('Posiciones/2posicion{0}.npy'.format(it))
       v2 = np.load('Velocidades/2velocidad{0}.npy'.format(it))
       xi = np.load('Posiciones/posicionion{0}.npy'.format(it))
       vi = np.load('Velocidades/velocidadion{0}.npy'.format(it))

       plt.figure('Espacio de fase')
       plt.plot(xi,vi,'.',markersize=0.5,c='b')
       plt.plot(x1,v1,'.',markersize=0.5,c='m')
       plt.plot(x2,v2,'.',markersize=0.5,c='c')
#       plt.title(u'Espacio de fase',fontsize = 18)
       plt.xticks(size='larger')
       plt.yticks(size='larger')
#       plt.xlabel('x',fontsize = 13)
#       plt.ylabel('v',fontsize = 13)
       plt.tick_params(labelsize=16)
       
#       plt.ylim(-5,5)
#       plt.xlim(0,120)
       plt.savefig('EspacioDeFase/EspacioDeFase{0}.png'.format(it))
       plt.show()
       plt.draw() 
       plt.clf()
       
       v=sp.append(v1,v2)
       OD.fdv(Npvel,v0,dv,ve,v,npe,L,it)

#
#GRAFICA LOS CAMPOS EN LOS TIEMPOS tfield
for i in range (0, tiempo):
    if i % tfield == 0:
        plt.figure('Potencial en el tiempo')
        plt.plot(X,phi[:,i],c='r')
#       plt.title(u'Perfil espacial Potencial Eléctrico',fontsize = 18)
        plt.xticks(size='larger')
        plt.yticks(size='larger')
#       plt.xlabel('x',fontsize = 13)
#       plt.ylabel('$\phi$',fontsize = 13)
        plt.xticks(size='larger')
        plt.yticks(size='larger')
        plt.tick_params(labelsize=16)
        plt.ylim(-6,6)
        #plt.xlim(0,1024)
        plt.savefig('Graficas/Potencialperfil{0}.png'.format(i))
        plt.figure('Campo en el tiempo')
        plt.plot(X,campo[:,i],c='r')
#       plt.title(u'Perfil espacial Potencial Eléctrico',fontsize = 18)
        plt.xticks(size='larger')
        plt.yticks(size='larger')
#       plt.xlabel('x',fontsize = 13)
#       plt.ylabel('$\phi$',fontsize = 13)
        plt.xticks(size='larger')
        plt.yticks(size='larger')
        plt.tick_params(labelsize=16)
        plt.ylim(-6,6)
        #plt.xlim(0,1024)
        plt.savefig('Graficas/Campoperfil{0}.png'.format(i))
        plt.figure('Densidad en el tiempo')
        plt.plot(X,rho[:,i],c='r')
#       plt.title(u'Perfil espacial Potencial Eléctrico',fontsize = 18)
        plt.xticks(size='larger')
        plt.yticks(size='larger')
#       plt.xlabel('x',fontsize = 13)
#       plt.ylabel('$\phi$',fontsize = 13)
        plt.xticks(size='larger')
        plt.yticks(size='larger')
        plt.tick_params(labelsize=16)
        plt.ylim(-6,6)
        #plt.xlim(0,1024)
        plt.savefig('Graficas/Densidadperfil{0}.png'.format(i))


#GRAFICAR LA ENERGIA EN EL TIEMPO
plt.figure('Energia en el tiempo')
plt.plot(T,EK,c='r')
plt.plot(T,EP,c='b')
plt.xticks(size='larger')
plt.yticks(size='larger')
#      plt.xlabel('x',fontsize = 13)
plt.xticks(size='larger')
plt.yticks(size='larger')
plt.tick_params(labelsize=16)
plt.ylim(-6,6)
#plt.xlim(0,1024)
plt.savefig('Graficas/Energiaperfil{0}.png'.format(tiempo))
plt.show()
plt.draw()
plt.clf()
     

#GRAFICA RELACION DE DISPERSION 
OD.relacion_dispersion(N,tiempo,Np,campo,dt,L)
