import numpy as np
import math


L2 = lambda v : np.sqrt(1./np.dot(*v.shape)*np.einsum('ij->', (np.abs(v))**2))
Linf = lambda v : np.max(np.abs(v))

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


def FromDat(x, y, Re, **kwargs):
    field = np.genfromtxt(kwargs.get('name', 1.))[1:,:]
    if(type(field)==np.ndarray):
        if(field.shape==(len(x), len(y))):
            return field
        else:
            print("Specified velocity field does not match grid initialized.")
            

def TaylorGreen(x, y, Re, **kwargs):
    kappa = kwargs.get('kappa', 1.)
    t = kwargs.get('time', 0.)
    field = 2 * kappa * np.cos(kappa * x) * np.cos(kappa * y[:, np.newaxis]) *\
                    np.exp(-2 * kappa**2 * t / Re)
    return field


def ShearLayer(x, y, Re, **kwargs):
    delta = kwargs.get('delta', 0.005)
    sigma = kwargs.get('sigma', 15./np.pi)
    field = delta * np.cos(x) - sigma * np.cosh(sigma * (y[:,np.newaxis] -\
                                                0.5*np.pi))**(-2)
    field+= delta * np.cos(x) + sigma * np.cosh(sigma * (1.5*np.pi -\
                                                y[:,np.newaxis]))**(-2)
    return field


def ConvectiveVortex(x, y, Re, **kwargs):
    Uinf = kwargs.get('Uinf', 1.)
    beta = kwargs.get('beta', 1./50.)
    R = kwargs.get('R', 0.005*np.pi)
    dx=x[1]-x[0]; dy=y[1]-y[0]
    # radial distance to vortex core
    rx = x - np.pi
    ry = y[:,np.newaxis]-np.pi
    r = np.sqrt(rx**2+ry**2)

    # init field
    # u = Uinf*(1-beta*(y[:,np.newaxis]-np.pi)/R*np.exp(-r**2/2))
    # v = Uinf*beta*(x-np.pi)/R*np.exp(-r**2/2)
    beta = 5.
    u = Uinf-beta/(2*np.pi)*np.exp(0.5*(1-r**2))*ry
    v = Uinf+beta/(2*np.pi)*np.exp(0.5*(1-r**2))*rx
    return Curl(u, v, dx, dy)


def McWilliams(x, y, Re, **kwargs):
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


def EnergySpectrum(k, s=3, kp=12):

    # normalise the spectrum
    a_s = (2*s + 1)**(s + 1) / (2**s * math.factorial(s))

    # compute sectrum at this wave number
    E = a_s/ (2 * kp) * (k / kp)**(2*s + 1) * np.exp(-(s + 0.5) * (k / kp)**2)

    return E


def PhaseFunction(kx, ky):

    # half the array size
    lenx2 = len(kx)//2; leny2 = len(ky)//2

    # define phase array
    xi = np.zeros((len(kx), len(ky)))

    # compute phase field in k space, need more points because of how 
    # python organises wavenumbers
    zeta = 2 * np.pi * np.random.rand(lenx2+1, leny2+1)
    eta  = 2 * np.pi * np.random.rand(lenx2+1, leny2+1)

    # quadrant \xi(kx,ky) = \zeta(kx,ky) + \eta(kx,ky)
    xi[:lenx2, :leny2] = zeta[:-1,:-1] + eta[:-1,:-1]

    # quadrant \xi(-kx,ky) = -\zeta(kx,ky) + \eta(kx,ky)
    xi[lenx2:, :leny2] = np.flip(-zeta[1:,:-1] + eta[1:,:-1], 0)

    # quadrant \xi(-kx,-ky) = -\zeta(kx,ky) - \eta(kx,ky)
    xi[lenx2:, leny2:] = np.flip(-zeta[1:,1:] - eta[1:,1:])

    # quadrant \xi(kx,-ky) = \zeta(kx,ky) - \eta(kx,ky)
    xi[:lenx2, leny2:] = np.flip(zeta[:-1,1:] - eta[:-1,1:], 1)

    return np.exp(1j * xi)


def DecayingTurbulence(x, y, Re, **kwargs):
    """
    Generates random vorticity field, see:
        Omer San, Anne E. Staples : High-order methods for decaying two-dimensional homogeneous isotropic turbulence
    """
    # Fourier mesh
    nx = len(x); kx = np.fft.fftfreq(nx, d=1./nx)
    ny = len(y); ky = np.fft.fftfreq(ny, d=1./ny)

    # define 2D spectrum array
    w_hat = np.empty((len(kx), len(ky)), dtype=np.complex128)

    # compute vorticity field in k space
    k = np.sqrt(kx**2 + ky[:, np.newaxis]**2)
    w_hat = np.sqrt((k / np.pi) * EnergySpectrum(k))

    # add random phase
    whatk = w_hat * PhaseFunction(kx, ky)

    # transforms initial field in physical space
    w = np.fft.ifft2(whatk) * nx * ny

    return np.real(w)