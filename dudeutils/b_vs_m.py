import numpy as np
H0=70.
omega_m=0.27
omega_b=0.044
omega_de=0.73
omega_k=0.0
sigma_8=0.9


def dX_(z,delz):#tytler 2009
    return delz*(1+z)/np.sqrt((1+z)*(1+z*omega_m)-z*(2+z)*omega_de)

#def dX(z,delz):#eqn from peroux et al 2003
#    return delz*(1+z)**2./np.sqrt((1.+z)**2. *(1+z*omega_m)-z*(2.+z)*omega_de)


def X(zf,delz):
    z=0.
    X=0.
    while z<zf:
        X+=dX(z,delz)
        z+=delz
    return X

def H(z):
    return H0*np.sqrt(omega_m*(1+z)**3.+omega_k*(1+z)**2. +omega_de)

def f(ablst, Nrng):
    _abs=[it for it in ablst if Nrng[0]<it.N<=Nrng[1]]
    num=float(len(_abs))
    redshift_rng=np.sum([it.X for it in _abs])
    delN=10.**Nrng[1] - 10.**Nrng[0]
    return num/ (delN*redshift_rng)
