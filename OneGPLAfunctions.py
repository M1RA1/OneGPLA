import math
import matplotlib.pyplot as plt
import scipy as sp
import numpy
from scipy.fftpack import fft2  # FFT dos dimensional
from scipy.fftpack import fftshift as shift  # Corrimiento al cero


def Maxwell(v,vd,vt):
    Distribucion = sp.exp(-(v - vd) * (v - vd) / (2* (vt*vt)))
    return Distribucion
def CargarRandomPosicion(Np,L):
    pst= sp.zeros(Np) #pst son las posiciones del sistema
    for i in range (Np):
        pst[i] = sp.random.uniform(0, L)
    return pst
def CargarVelocidadMaxwell1D(Np, n, vd, vt ):
    v = sp.zeros(Np)
    C = n / ((vt)*((2 * sp.pi)**(1 / 2))) #Constante de Normalizacion
    fmax = C
    vmin1 = vd - 5.0 * vt
    vmax1 = vd + 5.0 * vt
    for i in range(Np):
        while True:
             vtemp1 = vmin1 + (vmax1 - vmin1) * (sp.random.random())
             f = C* Maxwell(vtemp1, vd, vt)
             xtempl = fmax * (sp.random.random())
             if xtempl < f:
                  break
        v[i] = vtemp1
    return v
def CondicionesDeFrontera(Np,pst, L):
    #  Periodicos
    # pst son las posiciones
    for i in range(Np):
        if (pst[i] < 0.0):
            while(pst[i] < 0.0):
                  pst[i] += L
        elif (pst[i] >= L):
            while(pst[i] >=L):
                  pst[i] -= L

    return pst
def densidad(rho, dh, parts):
    np = len(parts)
    ni = len(rho)

    for i in range(0, ni):
                rho[i] = 0
    # Esparcir las partículas en la rejilla
    for p in range(0, np):
    # Adquirimos la posición de la partícula en la rejilla
             P = parts[p]
             x, v, tag, q, m, qm = P
             lc = x / dh
             i = math.floor(lc)
             d = lc - i
             rho[i] = rho[i] + q * (1 - d)
             rho[i + 1] = rho[i + 1] + q * d

    rho[ni - 1] = rho[ni - 1] + rho[0]
    rho[0] = rho[ni - 1]

     # Divide por el tamaño de las celdas para la densidad de carga
    for i in range(0, ni):
            rho[i] = rho[i] / dh

    # Remover ruido
    for i in range(0, ni):
            if (math.fabs(rho[i]) < 1e-10):
                rho[i] = 0
    return rho
def potentialGaussSeidel(ni, rho, phi, tol, dh, EPS0, it):
    phi[0]=0
    s=0
    for j in range(0,ni):
        phi[j]=0

    for solver_it in range(0,20000):
        for i in range(0,ni-1):
            im = i-1
            if im<0:
                im=ni-2
            ip=i+1
            if ip==ni-1:
                ip=0
            g = 0.5*((rho[i]/EPS0)*dh*dh+phi[im]+phi[ip])
            phi[i] = phi[i]+1.4*(g-phi[i])


        #Verifica convergencia
        if solver_it%25==0:
            sum=0
            for i in range(1,ni-1):
                ip = i+1
                if (ip==ni-1):
                    ip=0
                res=rho[i]/EPS0+(phi[i-1]-2*phi[i]+phi[ip])/(dh*dh)
                sum = sum+res*res

            l2=math.sqrt(sum/ni)/1
            if solver_it % 200 == 0:
                if it%5==0:
                   print("Solver", l2, tol)

            if l2<tol:
                s=1
                print('RUPTURA it', it, 'solver_it', solver_it)
                break;
    if s==0:
        print("NO CONVERGIO")
    phi[ni-1]=phi[0]
    return phi
def potentialFourier(Nx, dx, rho, EPS0):
    rho_k = numpy.fft.fftn(rho)
    Wx = numpy.exp(2 * 1j * numpy.pi / Nx)
    Wn = 1.0 + 0.0j
    dx_2 = dx * dx

    for n in range(Nx):
        denom = (2 - Wn - 1.0 / Wn)*EPS0
        if denom != 0:
            rho_k[n] *= dx_2 / denom
        Wn *= Wx

    phi = numpy.fft.ifftn(rho_k)
    phi = numpy.real(phi)
    return phi
def campo(ni, phi, ef, dh):
       for i in range(0, ni):
           im = i - 1
           ip = i + 1
           if im < 0:
               im = ni - 2
           if ip > ni - 1:
               ip = 1
           ef[i] = (phi[im] - phi[ip]) / (2 * dh)
       return ef
def movimiento( ef, dt, dh, parts, L, it):
    np = len(parts)
    for p in range(0, np):
        P = parts[p]
        x, v, tag, q, m, qm = P
        # Interpolar el campo a la posición de la partícula
        lc = x / dh
        i = math.floor(lc)
        d = lc - i
        ef_p = ef[i] * (1 - d) + ef[i + 1] * d
        #Primer paso, método leap frog
        if it==0:
           v=v-0.5*ef_p*qm*dt
        # Cargar velocidades
        v = v + ef_p * qm * dt
        # update position
        x = x + v * dt
        # Fronteras periódicas para las partículas
        if x < 0:
            while(x<0):
                  x = x + L
        if x >= L:
            while(x>=L):
                  x = x - L
        parts[p] = x, v, tag, q, m, qm  # Retorna los valores cambiados en STEP
    return parts
def energias(rho,phi,parts):
    np = len(parts)
    ni = len(rho)
    KE = 0
    PE = 0
    for p in range(0, np):
        P = parts[p]
        x, v, tag, q, m, qm = P
        KE = KE + m * v * v
    KE = KE * 0.5

    for i in range(0, ni):
        PE = PE + rho[i] * phi[i]
    PE = PE * 0.5
    return KE, PE
def MaxMin(var):
    # Se utiliza para cargar el valor  máximo y mínimo de alguna magnitud
    ni=len(var)
    var_min = var[0]
    var_max = var[0]
    for i in range(0, ni):
        if var[i] < var_min:
            var_min = var[i]
        if var[i] > var_max:
            var_max = var[i]
    print('min =', var_min)
    print('max =', var_max)
    return var_min, var_max

def stepion(ni, rho, rhoi, phi, ef, dt, dh, parts, partsi, EPS0, L, it):
    np = len(parts)
    rho = densidad(rho, dh, parts)

    rhoi = densidad(rhoi, dh, partsi)
    for i in range(0, ni):
        rho[i] = rho[i] + rhoi[i]
    # Potencial por medio de Fourier
    phi = potentialFourier(ni, dh, rho, EPS0)
    # Campo eléctrico
    ef = campo(ni, phi, ef, dh)
    #Cargar posiciones
    parts = movimiento(ef, dt, dh, parts, L, it)
    # partsi = PIC.movimiento(ef, dt, dh, partsi, L, it)
    return rho, phi, ef, parts
def stepmoveion(ni, rho, rhoi, phi, ef, dt, dh, parts, partsi, EPS0, L, it):
    np = len(parts)
    rho = densidad(rho, dh, parts)
    rhoi = densidad(rhoi, dh, partsi)
    for i in range(0, ni):
        rho[i] = rho[i] + rhoi[i]
    # Potencial por medio de Fourier
    phi = potentialFourier(ni, dh, rho, EPS0)
    # Campo eléctrico
    ef = campo(ni, phi, ef, dh)
    # Cargar posiciones
    partsi = movimiento(ef, dt, dh, partsi, L, it)
    parts = movimiento(ef, dt, dh, parts, L, it)
    return rho, phi, ef, parts, partsi

def plots(parts, L, vmag0):
    np = len(parts)
    plt.xlim(0, L)
    plt.ylim(-3 * vmag0, 3 * vmag0)
    EXB = []
    EVB = []
    EXR = []
    EVR = []
    for p in range(0, np):
        P = parts[p]
        x, v, tag, q, m, qm = P
        if tag == 0:
            EXR.append(x)
            EVR.append(v)
        else:
            EXB.append(x)
            EVB.append(v)

            # dibuja gráfico rojo
    plt.scatter(EXR, EVR)
    # dibuja gráfico azul
    plt.scatter(EXB, EVB)
    plt.pause(0.001)
    plt.clf()
    # plt.show()
def nextpow2(longitud_malla):
    """Potencia más cercana de dos"""
    n = 1
    while n < longitud_malla: n *= 2
    return n

def fdv(Npvel,v0,dv,ve,v,npe,L,it):
    FD = sp.zeros(Npvel)
    for i in range(npe):
    #get particle logical coordinate*/
#        P=parts[i]
#        x,v,tag,q,m,qm=P 
        try:
           g = (v[i] + v0) / dv
           jp = int(g)
           f1 = g - jp
           f2 = 1 - f1
           FD[jp] = FD[jp] + f2
           FD[jp + 1] = FD[jp + 1] + f1
        except:
            print('Particula se escapó :c')
        
    FD=FD/L
    plt.figure('Función de Distribución')
    plt.plot(ve, FD,c='darkmagenta')
#    plt.title(u'Espacio de velocidades',fontsize = 18)
    plt.xticks(size='larger')
    plt.yticks(size='larger')
#    plt.xlabel('v',fontsize = 13)
#    plt.ylabel('f(v)',fontsize = 13)
    plt.tick_params(labelsize=16)

#    plt.ylim(0,10000)
#    plt.xlim(-5,5)
    plt.savefig('FuncionDistribucion/FuncionDistribucion{0}.png'.format(it))
    plt.pause(0.001)
    plt.clf()
def relacion_dispersion(npuntos_malla,npasos_temporales,nparticulas,E_acumulador,dt,longitud_malla):

    # Se calcula la frecuencia angular y espacial minima y maxima (Ver codigo KEMPO1):
    omega_min = 2*sp.pi/(dt)/2/(npasos_temporales / 2)
    omega_max = omega_min * (npasos_temporales / 2)
    k_min = 2 * sp.pi / (npuntos_malla)
    k_max = k_min * ((npuntos_malla / 2)-1)
    # Se crean los vectores de frecuencias espacial y angular simuladas:
    nxtp2=nextpow2(npuntos_malla)
    k_simulada = sp.linspace(-k_max,k_max,nxtp2)*10
    omega_simulada = sp.linspace(-omega_max,omega_max,nxtp2)/10 
    #El diez es normalizando con respecto a vt
    # Se genera una matriz de frecuencias angular y espacial:
    K, W = sp.meshgrid(k_simulada,omega_simulada)     
    # Se muestrea la matriz espacio-temporal:
    E_acumulador_muestreado = E_acumulador[0:npuntos_malla:1,0:npasos_temporales:5]    
    # Se efectua la FFT sobre la matriz espacio temporal muestreada, luego el 
    # valor absoluto de dicha matriz y el corrimiento al cero de las frecuencias:
    E_wk = fft2(E_acumulador_muestreado,(nextpow2(npuntos_malla),nextpow2(npuntos_malla)))/longitud_malla
    E_wk_absoluto = abs(E_wk)
    E_wk_shift = shift(E_wk_absoluto)    
#    plt.xticks(np.linspace(0,0.7,6), fontsize = 18)
#    plt.yticks(np.linspace(0,1,5), fontsize = 18)
    # Se grafica la relacion de dispersion simulada:
    plt.contourf(-K,W,E_wk_shift, 8, alpha=.75, cmap='rainbow')
#    plt.xlim(0,0.7)
#    plt.ylim(0.0,1.1)
    clb=plt.colorbar()
    clb.ax.set_title('|E|')
    plt.savefig('Graficas/Relaciondispersion.png')