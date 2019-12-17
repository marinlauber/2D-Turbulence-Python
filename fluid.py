#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Marin Lauber"
__copyright__ = "Copyright 2019, Marin Lauber"
__license__ = "GPL"
__version__ = "1.0.1"
__email__  = "M.Lauber@soton.ac.uk"

import numpy as np
import matplotlib.pyplot as plt
# import pyfftw

class Fluid(object):

    def __init__(self, nx, ny, Re, dt=0.0001, pad=3./2.):

        # input data
        self.nx = nx
        self.ny = ny; self.nk = self.ny//2+1
        self.Re = Re; self.ReI = 1./self.Re
        self.dt = dt
        self.pad = pad
        self.time = 0.
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
        self.w0 = np.empty((self.nx,self.nk))
        self.dwdt = np.empty((self.nx,self.ny))

    
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
        KEaux = self._spec_variance(psih)
        psi = psih/np.sqrt(KEaux)

        # inverse Laplacian in k-space
        wh = self.k2 * psi
        
        # vorticity in physical space
        self.w = np.fft.irfft2(wh)
        self.w = self.w/(0.5*abs(self.w).max())-1.


    def init_solver(self):

        try:    
            self.k2
        except AttributeError:
            self.k2 = self.kx[:self.nk]**2 + self.ky[:,np.newaxis]**2
            self.fk = self.k2 != 0.0

        # initialise array required for solving
        self.wh = np.zeros((self.nx,self.nk), dtype=np.complex128)
        self.psih = np.zeros((self.nx,self.nk), dtype=np.complex128)
        self.dwhdt = np.zeros((self.nx,self.nk), dtype=np.complex128)

        # utils
        self.mx = int(self.pad * self.nx)
        self.my = int(self.pad * self.nk)
        self.padder = np.ones(self.mx, dtype=bool)
        self.padder[int(self.nx/2):int(self.nx*(self.pad-0.5)):] = False

        # populate those arrays
        self.w0 = self.w

        # ṣpectral filter
        try:
            self.fltr
        except AttributeError:
            self._init_filter()


    def w_to_wh(self):
        self.wh = np.fft.rfft2(self.w, axes=(-2,-1))
        

    def wh_to_w(self):
        self.w = np.fft.irfft2(self.wh, axes=(-2,-1))


    def dwhdt_to_dwdt(self):
        self.dwdt = np.fft.irfft2(self.dwhdt, axes=(-2,-1))
        

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
        self.u = np.fft.irfft2(self.ky[:,np.newaxis]*self.psih)


    def get_v(self):
        """
        Spectral differentiation to get:
            v = -d/dx \psi
        """
        self.v = -np.fft.irfft2(self.kx[:self.nk]*self.psih)


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
        # self.dt = 0.005 * np.min(np.hstack((self.dx, self.dy))) /\
        #  np.max(np.hstack((self.u, self.v)))
        # self.dt = np.minimum(self.dt, 0.0001)


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
            self._get_psih()

            # get convective forces (resets dwhdt)
            self._add_convection()

            # add diffusion (changes to C-space)
            self._add_diffusion()

            # step in time
            self.w = self.w0 + (self.dt/k) * self.dwdt

        self.time += self.dt
        self._cfl_limit()
    

    def _get_psih(self):
        """
        Spectral stream-function from spectral vorticity
            hat{\psi} = \hat{\omega} / (k_x^2 + k_y^2)
        """
        self.w_to_wh()
        self.psih[self.fk] = self.wh[self.fk] / self.k2[self.fk]


    def _add_diffusion(self):
        """
        Diffusion term of the Navier-Stokes
            D = p * 1/Re * (-k_x^2 -k_y^2) * \hat{\omega}
        
        Note: This resets the value in self.w when called
              The penalty value is required for the RK method
        """
        self.dwhdt -= self.ReI*self.k2*self.wh
        self.dwhdt_to_dwdt()


    def _add_convection(self):
        """
        Convective term
            N = -d/dy \psi * d/dx \omega + d/dx \psi * d/dy \omega
        To prevent alliasing, we zero-pad the array before using the
        convolution theorem to evaluate it in physical space.
        """
        self.dwhdt = -1*self._convolve(1j*self.ky[:,np.newaxis]*self.psih,
                                       1j*self.kx[:self.nk]*self.wh)
        self.dwhdt +=   self._convolve(1j*self.kx[:self.nk]*self.psih,
                                       1j*self.ky[:,np.newaxis]*self.wh)
        # self._add_spec_filter()


    def _convolve(self, a, b):
        """
        Evaluate convolution sum. This involves three transforms
        """
        # zero-padded temp arrays
        tmp = np.zeros((self.mx,self.my), dtype=np.complex128)
        tmp[self.padder, :self.nk] = a
    
        # fft with these new coeff, padded with zeros
        r = np.fft.irfft2(tmp, axes=(-2,-1))
        tmp *= 0.0j; tmp[self.padder, :self.nk] = b

        # multiplication in physical space, this saves one temp array
        r *= np.fft.irfft2(tmp, axes=(-2,-1))*self.pad**(2)
        tmp = np.fft.rfft2(r, axes=(-2,-1))
        
        return tmp[self.padder, :self.nk] # truncate fourier modes
    

    def _add_spec_filter(self):
        self.dwhdt *= self.fltr

    
    def _spec_variance(self, ph):
        # spectral filter
        self._init_filter()

        # only half the spectrum for real ffts, needs spectral normalisation
        var_dens = 2 * np.abs(self.fltr*np.sqrt(self.k2)*ph)**2 / (self.nx*self.ny)**2
        var_dens[..., 0] /= 2
        var_dens[...,-1] /= 2

        return var_dens.sum(axis=(-2,-1))


    def _tke(self):
        ke = .5*self._spec_variance(np.sqrt(self.k2)*self.psih)
        return ke.sum()


    def _compute_spectrum(self, res):
        self._get_psih()
        # angle averaged TKE spectrum
        tke = np.real(.5*self.k2*self.psih*np.conj(self.psih))
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
        np.savetxt(str(folder)+"vort_"+str("%06d"%iter)+".dat", np.vstack((s, self.w)))


    def show_vort(self):
        p=plt.imshow(self.w, cmap="RdBu")
        plt.colorbar(p)
        plt.xticks([]); plt.yticks([])
        plt.show()
    

    def show_spec(self):
        plt.figure()
        plt.imshow(np.real(self.wh))
        plt.show()


    # def show_vel(self):
    #     if(self.uptodate!=True):
    #         self.w_to_wh()
    #         self._get_psih()
    #         self.get_u()
    #         self.get_v()
    #     plt.figure()
    #     plt.quiver(self.x, self.y, self.u, self.v)
    #     plt.xlabel("x"); plt.ylabel("y")
    #     plt.show()


    def run_live(self, stop, every=100):
        from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
        iterr = 0
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
            iterr += 1
            if(iterr % every == 0):
                im.set_data(self.w)
                fig.canvas.draw()
                fig.canvas.flush_events()
                print("Iteration \t %d, time \t %f, time remaining \t %f. TKE: %f" %(iterr,
                      self.time, stop-self.time, self._tke()))

 

# if __name__=="__main__":
#     flow = Fluid(128, 128, 1)
#     flow.init_field("TG")
#     flow.init_solver()
#     print(flow._tke())
