import numpy as np


def _spec_variance(ph):
        # only half the spectrum for real ffts, needs spectral normalisation
        nx, nk = ph.shape
        ny = (nk-1)*2
        var_dens = 2 * np.abs(ph)**2 / (nx*ny)**2
        # only half of coefs [0] and [nx/2+1] due to symmetry in real fft2
        var_dens[..., 0] /= 2.
        var_dens[...,-1] /= 2.

        return var_dens.sum(axis=(-2,-1))


def Curl(u, v, dx, dy):
    curl = np.zeros_like(u)
    curl[1:-1,1:-1] = (u[1:-1,2:]-u[1:-1,:-2])/dy - ((v[2:,1:-1]-v[:-2,1:-1])/dx)
    return curl


def TaylorGreen(x, y, Re=1., t=0.0, kappa=2.):
    field = 2 * kappa * np.cos(kappa * x) * np.cos(kappa * y[:, np.newaxis]) *\
                    np.exp(-2 * kappa**2 * t / Re)
    return field


def ShearLayer(x, y, delta=0.005, sigma= 15./np.pi):
    field = delta * np.cos(x) - sigma * np.cosh(sigma * (y[:,np.newaxis] -\
                                                0.5*np.pi))**(-2)
    field+= delta * np.cos(x) + sigma * np.cosh(sigma * (1.5*np.pi -\
                                                y[:,np.newaxis]))**(-2)
    return field


def ConvectiveVortex(x, y, Uinf=1., beta=1./50., R=0.005*np.pi):
    dx=x[1]-x[0]; dy=y[1]-y[0]
    # radial distamce to vortex core
    r = np.sqrt((x-np.pi)**2+(y[:,np.newaxis]-np.pi)**2)

    # init field
    u = Uinf*(1-beta*(y[:,np.newaxis]-np.pi)/R*np.exp(-r**2/2))
    v = Uinf*beta*(x-np.pi)/R*np.exp(-r**2/2)

    return Curl(u, v, dx, dy)


def McWilliams(x, y):
    """
    Generates McWilliams vorticity field, see:
        McWilliams (1984), "The emergence of isolated coherent vortices in turbulent flow"
    """
    # Fourier mesh
    nx = len(x); kx = np.fft.fftfreq(nx, d=1./nx)
    ny = len(y); ky = np.fft.fftfreq(ny, d=1./ny)
    nk = ny//2+1

    # generate variable
    k2 = kx[:nk]**2 + ky[:,np.newaxis]**2
    fk = k2 != 0.0

    # ensemble variance proportional to the prescribed scalar wavenumber function
    ck = np.zeros((nx, nk))
    ck[fk] = (np.sqrt(k2[fk])*(1+(k2[fk]/36)**2))**(-1)
    
    # Gaussian random realization for each of the Fourier components of psi
    psih = np.random.randn(nx, nk)*ck+\
            1j*np.random.randn(nx, nk)*ck

    # ṃake sure the stream function has zero mean
    cphi = 0.65*np.max(kx)
    wvx = np.sqrt(k2)
    filtr = np.exp(-23.6*(wvx-cphi)**4.)
    filtr[wvx<=cphi] = 1.
    KEaux = _spec_variance(filtr*np.sqrt(k2)*psih)
    psi = psih/np.sqrt(KEaux)

    # inverse Laplacian in k-space
    wh = k2 * psi
    
    # vorticity in physical space
    field = np.fft.irfft2(wh)
    return field
