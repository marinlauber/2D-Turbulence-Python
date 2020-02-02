#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Marin Lauber"
__copyright__ = "Copyright 2019, Marin Lauber"
__license__ = "GPL"
__version__ = "1.0.1"
__email__  = "M.Lauber@soton.ac.uk"

import numpy as np
import matplotlib.pyplot as plt

class Fluid(object):

    def __init__(self, nx, ny, Re, dt=0.0001):

        # input data
        self.nx = nx
        self.ny = ny; self.nk = self.ny//2+1
        self.Re = Re; self.ReI = 1./self.Re
        self.dt = dt
        self.time = 0.
        self.it =0
        self.uptodate = False
        self.filterfac = 23.6

        # we assume 2pi periodic domain in each dimensions
        self.x = np.linspace(0, 2*np.pi, nx, endpoint=False)
        self.y = np.linspace(0, 2*np.pi, ny, endpoint=False)

        # physical grid
        self.dx = 2*np.pi/self.nx
        self.dy = 2*np.pi/self.ny

        # fourier grid
        self.kx = np.fft.fftfreq(self.nx)*self.nx
        self.ky = np.fft.fftfreq(self.ny)*self.ny
        
        # different attribute
        self.u = np.empty((self.nx,self.ny))
        self.v = np.empty((self.nx,self.ny))
        self.w = np.empty((self.nx,self.ny))

    
    def init_field(self, field="Taylor-Green", t=0.0, kappa=2., delta=0.005, sigma= 15./np.pi):
        """
        Inital flow field following the Taylor-Green solution of the Navier-Stokes
        or a double shear layer.

            Params:
                field: string
                    -Type of field to initialise
        """
        if(type(field)==str):
            if(field=="TG" or field=="Taylor-Green"):
                self.w = 2 * kappa * np.cos(kappa * self.x) * np.cos(kappa * self.y[:, np.newaxis]) *\
                        np.exp(-2 * kappa**2 * t / self.Re)
            elif(field=="SL" or field=="Shear Layer"):
                self.w = delta * np.cos(self.x) - sigma * np.cosh(sigma * (self.y[:,np.newaxis] -\
                        0.5*np.pi))**(-2)
                self.w += delta * np.cos(self.x) + sigma * np.cosh(sigma * (1.5*np.pi -\
                        self.y[:,np.newaxis]))**(-2)
            elif(field=="McWilliams" or field=="MW84"):
                self.McWilliams1984()
            else:
                print("The specified field type %s is unknown.\nAvailable initial fields are"+\
                      ": \"Taylor-Green\", \"Shear Layer\"." % field)
        elif(type(field)==np.ndarray):
            if(field.shape==(self.nx,self.ny)):
                self.w = field
            else:
                print("Specified velocity field does not match grid initialized.")

        # populate array
        self._get_psi()
        self.w0 = self.w


    def init_solver(self, scheme="PADE-6"):

         # initialise array required for solving
        self.w0 = np.empty((self.nx,self.ny), dtype=np.float64)
        self.psi = np.zeros((self.nx,self.ny), dtype=np.float64)
        self.dwdt = np.zeros((self.nx,self.ny), dtype=np.float64)

        # pade schemes coefficients
        if scheme=="CDS-2":
            self.alpha = 0.; self.beta  = 1./(2.*self.dx); self.gamma = 0.
        if scheme == "CDS-4":
            self.alpha = 0.; self.beta  = 2./(3.*self.dx); self.gamma = -1./(12.*self.dx)
        if scheme == "PADE-4":
            self.alpha = 0.25; self.beta  = 3./(4.*self.dx); self.gamma = 0.
        if scheme == "PADE-6":
            self.alpha = 1./3.; self.beta  = 14./(18.*self.dx); self.gamma = 1./(36.*self.dx)
        
        # utils
        self.k2 = self.kx[:self.nk]**2 + self.ky[:,np.newaxis]**2
        self.fk = self.k2 != 0.0


    def McWilliams1984(self):

        # generate variable
        self.k2 = self.kx[:self.nk]**2 + self.ky[:,np.newaxis]**2
        self.fk = self.k2 != 0.0

        # emsemble variance proportional to the prescribed scalar wavenumber function
        ck = np.zeros((self.nx, self.nk))
        ck[self.fk] = np.sqrt(self.k2[self.fk]*(1+(self.k2[self.fk]/36)**2))**(-1)
        
        # Gaussian random realization for each of the Fourier components of psi
        psih = np.random.randn(self.nx, self.nk)*ck+\
               1j*np.random.randn(self.nx, self.nk)*ck

        # ṃake sure the stream function has zero mean
        psi = np.fft.irfft2(psih)
        psih = np.fft.rfft2(psi-psi.mean())
        KEaux = self._spec_variance(self.fltr*np.sqrt(self.k2)*psih)
        psi = psih/np.sqrt(KEaux)

        # inverse Laplacian in k-space
        wh = self.k2 * psi
        
        # vorticity in physical space
        self.w = np.fft.irfft2(wh)


    def _init_filter(self):
        cphi = 0.65*np.max(self.kx)
        wvx = np.sqrt(self.k2)
        filtr = np.exp(-self.filterfac*(wvx-cphi)**4.)
        filtr[wvx<=cphi] = 1.
        self.fltr = filtr


    def get_u(self):
        """
        Spectral differentiation to get:
            u = d/dy \psi
        """
        self.u = np.fft.irfft2(self.ky[:,np.newaxis]*np.fft.rfft2(self.psi))


    def get_v(self):
        """
        Spectral differentiation to get:
            v = -d/dx \psi
        """
        self.v = -np.fft.irfft2(self.kx[:self.nk]*np.fft.rfft2(self.psi))


    def _cfl_limit(self):
        """
        Adjust time-step based on the courant condition
        
        Note: this assumes that you initial velocity field is correctly normalized.
        """
        self.get_u()
        self.get_v()
        Dc = np.max(np.pi*((1.+abs(self.u))/self.dx + (1.+abs(self.v))/self.dy))
        Dmu = np.max(np.pi**2*(self.dx**(-2) + self.dy**(-2)))
        self.dt = np.sqrt(3.) / (Dc + Dmu)


    def update(self, s=3):
        """
        Low-storage S-order Runge-Kutta method from Jameson, Schmidt and Turkel (1981)
        Input:
            s : float
                desired order of the method, default is 3rd order
        """
        # iniitalise field
        self.w0 = self.w

        for k in range(s, 0, -1):
            # invert Poisson equation for the stream function (changes to k-space)
            self._get_psi()

            # get convective forces (resets dwdt)
            self._add_convection()

            # add diffusion (changes to C-space)
            self._add_diffusion()

            # step in time
            self.w = self.w0 + (self.dt/k) * self.dwdt

        self.time += self.dt
        self._cfl_limit()
        self.it += 1
    

    def _get_psi(self):
        """
        Spectral stream-function from spectral vorticity
            hat{\psi} = \hat{\omega} / (k_x^2 + k_y^2)
        """
        wh = np.fft.rfft2(self.w, axes=(-2, -1))
        psih = np.zeros_like(wh)
        psih[self.fk] = wh[self.fk] / self.k2[self.fk]
        self.psi = np.fft.irfft2(psih, axes=(-2, -1))


    def _ddx(self, u, dx): return self._pade(u, dx)
    

    def _ddy(self, u, dy): return self._pade(u.T, dy).T


    def _pade(self, a, dx):
        """
        Second, fourth- and sixth-order compact scheme for first derivative of periodic field a:
        
        Paramaters:
            a     : array of float
                    values of the fields at the grid points the array is assumed periodic
                    such that the periodic point is NOT included
            dx    : float
                    grid spacing, must be constant
        Output:
            diff : array of float
                   first derivative of the field
        """
    
        # temp array
        b = np.empty_like(a)
        
        # compact scheme on interior points
        b[:, 2:-2] = self.beta * (a[:, 3:-1] - a[:, 1:-3]) + self.gamma * (a[:, 4:] - a[:, :-4])
        
        # boundary points
        b[:, -2] = self.beta * (a[:, -1] - a[:, -3]) + self.gamma * (a[: , 0] - a[:, -4])
        b[:, -1] = self.beta * (a[:,  0] - a[:, -2]) + self.gamma * (a[:,  1] - a[:, -3])
        b[:,  0] = self.beta * (a[:,  1] - a[:, -1]) + self.gamma * (a[:,  2] - a[:, -2])
        b[:,  1] = self.beta * (a[:,  2] - a[:,  0]) + self.gamma * (a[:,  3] - a[:, -1])
        
        # return
        if self.alpha == 0.:
            diff = b
        else:
            # build first row of circulant Pade coefficients matrix
            a = np.hstack((1., self.alpha, [0.*j for j in range(len(a)-3)], self.alpha))
            # solve tri-diagonal system using FFTs and circulant matrix properties
            diff = np.fft.irfft2(np.fft.rfft2(b)/np.fft.rfft(a))
        return diff


    def _add_convection(self):
        """
        Convective term
            -d/dy \psi * d/dx \omega + d/dx \psi * d/dy \omega
        To prevent alliasing, we zero-pad the array before using the
        convolution theorem to evaluate it in physical space.
        """
        self.dwdt = -1*self._ddy(self.psi, self.dy) * self._ddx(self.w, self.dx)
        self.dwdt +=   self._ddx(self.psi, self.dx) * self._ddy(self.w, self.dy)


    def _add_diffusion(self):
        """
        Diffusion term of the Navier-Stokes
            1/Re * (-k_x^2 -k_y^2) * \hat{\omega}
        """
        self.dwdt += self.ReI * ( self._ddx(self._ddx(self.w, self.dx), self.dx) +\
                                  self._ddy(self._ddy(self.w, self.dy), self.dy) )

    
    def _spec_variance(self, ph):
        # only half the spectrum for real ffts, needs spectral normalisation
        var_dens = 2 * np.abs(ph)**2 / (self.nx*self.ny)**2
        var_dens[..., 0] /= 2
        var_dens[...,-1] /= 2

        return var_dens.sum(axis=(-2,-1))


    def tke(self):
        psih = np.fft.rfft2(self.psi, axes=(-2,-1))
        ke = .5*self._spec_variance(np.sqrt(self.k2)*psih)
        return ke.sum()


    def enstrophy(self):
        eps = 0.5*abs(self.w)**2
        return eps.sum(axis=(-2,-1))


    def _compute_spectrum(self, res):
        psih = np.fft.rfft2(self.psi, axes=(-2,-1))
        # angle averaged TKE spectrum
        tke = np.real(.5*self.k2*psih*np.conj(psih))
        kmod = np.sqrt(self.k2)
        self.k = np.arange(1, self.nk, 1, dtype=np.float64) # niquist limit for this grid
        self.E = np.zeros_like(self.k)
        dk = (np.max(self.k)-np.min(self.k))/res

        #  binning energies with wavenumber modulus in threshold
        for i in range(len(self.k)):
            self.E[i] += np.sum(tke[(kmod<self.k[i]+dk) & (kmod>=self.k[i]-dk)])

    
    def plot_spec(self,res=200):
        self._compute_spectrum(200)
        plt.figure(figsize=(6,6))
        plt.plot(self.k, self.E, '-k', label="E(k)")
        plt.xlabel("k")
        plt.ylabel("E(k)")
        plt.legend()
        plt.show()


    def save_vort(self, folder, iter):
        s = np.zeros(self.ny); s[0]=self.time; s[1]=self.dt
        s[2] =self.tke(); s[3] = self.enstrophy()
        np.savetxt(str(folder)+"vort_"+str("%06d"%iter)+".dat", np.vstack((s, self.w)))


    def show_vort(self):
        p=plt.imshow(self.w, cmap="RdBu")
        plt.colorbar(p)
        plt.xticks([]); plt.yticks([])
        plt.show()
    

    def show_spec(self):
        plt.figure()
        plt.imshow(np.real(np.fft.fft2(self.w, axes=(-2,-1))))
        plt.show()


    def show_vel(self):
        self.get_u()
        self.get_v()
        plt.figure()
        plt.quiver(self.x, self.y, self.u, self.v)
        plt.xlabel("x"); plt.ylabel("y")
        plt.show()


    def run_live(self, stop, every=100):
        from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
        plt.ion()
        fig = plt.figure()
        ax = fig.add_subplot(111)
        im = ax.imshow(self.w, norm=None, cmap="RdBu")
        cax = make_axes_locatable(ax).append_axes("right", size="5%", pad="2%")
        cb = fig.colorbar(im, cax=cax)
        ax.set_xticks([]); ax.set_yticks([])
        while(self.time<=stop):
            #  update using RK
            self.update()
            if(self.it % every == 0):
                im.set_data(self.w)
                fig.canvas.draw()
                fig.canvas.flush_events()
                print(f"Iteration \t %d, time \t %f, time remaining \t %f. TKE: %f. ENS: %f" %(self.it,
                      self.time, stop-self.time, self.tke(), self.enstrophy()))


# if __name__=="__main__":
#     flow = Fluid(128, 128, 1)
#     flow.init_field("TG")
#     flow.init_solver()
#     print(flow._tke())
