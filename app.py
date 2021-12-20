import streamlit as st
import numpy as np
import scipy.integrate
import scipy.special
import scipy
import matplotlib.pyplot as plt


import numpy as np
import scipy.integrate
import scipy.special
import scipy
import matplotlib.pyplot as plt

#%%
bobth=[3.015e-3, 5.46e-3, 2.94e-3, 900, 1.32e-3, 6.03e-3]

args=[1e3,bobth,16.5e6,1,1000]

bo=args[1]

#%%
def cquad(func, a, b, **kwargs):
    """ funcion auxiliar"""
    def real_func(x):
        return np.real(func(x))
    def imag_func(x):
        return np.imag(func(x))
    real_integral = scipy.integrate.quad(real_func, a, b, **kwargs)
    imag_integral = scipy.integrate.quad(imag_func, a, b, **kwargs)
    return (real_integral[0] + 1j*imag_integral[0])

#%%
def ji(k,r1,r2):
    """ funcion auxiliar"""
    return scipy.integrate.quad(lambda x: x*scipy.special.jv(1,x) ,k*r1,k*r2)[0]

#%%
def expz(k,z1,z2):
    """ funcion auxiliar"""
    return ((np.exp(-k*z1)-np.exp(-k*z2)))/(k**(3))

#%%
def sigj(k,sigma,f,mur,z):
    """ funcion auxiliar"""
    mu0 = 4*np.pi*1e-7
    k2=1j*2*np.pi*f*mu0*mur*sigma
    la=np.sqrt(k**2+k2)
    return (2*k*mur*np.exp(la*z)/(la+k*mur))

#%%
def jhf(r, z, I, args):
    """ Calculo de densidad corriente sobre una placa semi-infinita """

    f=args[0]
    bo=args[1]
    sigma=args[2]
    mur=args[3]
    lmax=args[4]
    # jhf=list()
    mu0=4*np.pi*1e-7;
    r1=bo[0];
    r2=bo[1];
    dh=bo[2];
    N=bo[3];
    z1=bo[4];
    l0=bo[5]
    i0=N*I/((r2-r1)*dh)
    aint=1j*2*np.pi*f*sigma*mu0*i0
    inte=cquad(lambda k: scipy.special.j1(k*r)*ji(k,r1,r2)*expz(k,z1,z1+dh)*sigj(k,sigma,f,mur,z),0,lmax)
    return(aint*inte)

#%% Cálculo los vectores para graficar la dendidad de corriente en el primer
#   cuadrante

xv = np.arange(0, 0.011, 0.001)
yv = np.arange(0, 0.011, 0.001)
xmesh, ymesh = np.meshgrid(xv, yv) 

# Para graficar la corriente tengo que formar las matrices umesh y vmesh de 
# dimensión len(x)*len(y); que son las componentes de J


z = st.sidebar.slider('Densidad de corriente en el material a una distancia de z mm  de la superficie?', 0, 10, 0)


I = st.sidebar.slider('Corriente [mA]', 0.00, 0.10, 0.01)

Jt = np.zeros([len(xv), len(yv)])
Jx = np.zeros([len(xv), len(yv)])
Jy = np.zeros([len(xv), len(yv)])

for i,x in enumerate(xv): # Son las filas de las matrices Jr, Jrx y Jry
    for j,y in enumerate(yv): # Son la columnas de las matrices Jr, Jrx y Jry
        r = np.sqrt(x**2+y**2)
        
        
        # r,z,I, args
        jt = np.real(jhf(r,-z*1e-3,I,args))
        Jt[j,i] = jt
        
        if x != 0:
            t = np.arctan(y/x)
        else:
            t = np.pi / 2
        
        Jx[j,i] = -jt * np.sin(t)
        Jy[j,i] = jt * np.cos(t)


fig, (ax1, ax2) = plt.subplots(1, 2)
fig.suptitle('Horizontally stacked subplots')


ax1.quiver(xmesh, ymesh, Jx, Jy,scale=500000,units='inches')
ax1.quiver(xmesh, -ymesh, -Jx, Jy,scale=500000,units='inches')       
ax1.quiver(-xmesh, ymesh, Jx, -Jy,scale=500000,units='inches')    
ax1.quiver(-xmesh, -ymesh, -Jx, -Jy,scale=500000,units='inches')       
ax1.ticklabel_format(style='sci', scilimits=(-3,-3))







Jm=np.sqrt(Jx**2+ Jy**2)

@st.cache
def getmax():

    Jmax=Jm.max()

    return(Jmax)



Jm=np.sqrt(Jx**2+ Jy**2)/getmax()
colors='coolwarm'
#%% Grafica las lineas de contorno de la densidad de corriente 


cs0=ax2.contourf(xmesh,ymesh,Jm,
                  cmap = colors)
cs1=ax2.contourf(xmesh,-ymesh,Jm,
                  cmap = colors)
cs2=ax2.contourf(-xmesh,ymesh,Jm,
                  cmap = colors)
cs3=ax2.contourf(-xmesh,-ymesh,Jm,
                  cmap = colors)

plt.clabel(cs1, fmt='%2.1f', colors='k', fontsize=10,inline=True)
plt.ticklabel_format(style='sci', scilimits=(-3,-3))


fig.set_figheight(5)
fig.set_figwidth(10)

st.pyplot(fig) 
