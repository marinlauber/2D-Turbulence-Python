#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Marin Lauber"
__copyright__ = "Copyright 2019, Marin Lauber"
__license__ = "GPL"
__version__ = "1.0.1"
__email__  = "M.Lauber@soton.ac.uk"

import numpy as np
import matplotlib.pyplot as plt
import pyfftw

class Fluid(object):

    def __init__(self, nx, ny, Re, dt=0.0001, pad=3./2.):
        """
        initalizes the fluid, given a number or grid points in x and y. Sets flow parameters.
        Parameters:
            nx : intger
                - gird points in the x-direction
            ny : integer
                - grid points in the y-direction
            Re : float
                - Reynolds number of the flow
            dt : float
                - time-step, outdated as we use adaptive time-step
            pad : float
                - padding length for Jacobian evaluation
        """
        # input data
        self.nx = nx
        self.ny = ny; self.nk = self.ny//2+1
        self.Re = Re; self.ReI = 1./self.Re
        self.dt = dt
        self.pad = pad
        self.time = 0.
        self.uptodate = False
        self.filterfac = 23.6

        self.FFTW = True
        self.fftw_num_threads = 6

        # we assume 2pi periodic domain in each dimensions
        self.x, self.dx = np.linspace(0, 2*np.pi, nx, endpoint=False, retstep=True)
        self.y, self.dy = np.linspace(0, 2*np.pi, ny, endpoint=False, retstep=True)

        # fourier grid
        self.kx = np.fft.fftfreq(self.nx)*self.nx
        self.ky = np.fft.fftfreq(self.ny)*self.ny


    def init_solver(self):
        """
        Initalizes storage arrays and FFT objects. This is relatively expansive
        as pyFFTW find the fastest way to perform the transfoms, but it is only
        called once.
        """
        try:    
            self.k2
        except AttributeError:
            self.k2 = self.kx[:self.nk]**2 + self.ky[:,np.newaxis]**2
            self.fk = self.k2 != 0.0

        # utils
        self.mx = int(self.pad * self.nx)
        self.mk = int(self.pad * self.nk)
        self.my = int(self.pad * self.ny)

        # for easier slicing when padding
        self.padder = np.ones(self.mx, dtype=bool)
        self.padder[int(self.nx/2):int(self.nx*(self.pad-0.5)):] = False

        # initialise array required for solving
        self.u = self._empty_real()
        self.v = self._empty_real()
        self.w = self._empty_real()
        self.w0 = self._empty_real()
        self.dwdt = self._empty_real()

        self.uh = self._empty_imag()
        self.vh = self._empty_imag()
        self.wh = self._empty_imag()
        self.wh = self._empty_imag()
        self.psih = self._empty_imag()
        self.dwhdt = self._empty_imag()

        # assign padded arrays for non-linear term
        self.a = self._empty_imag((self.mx,self.mk))
        self.a1 = self._empty_imag((self.mx,self.mk))
        self.a2 = self._empty_imag((self.mx,self.mk))
        self.a3 = self._empty_imag((self.mx,self.mk))
        self.a4 = self._empty_imag((self.mx,self.mk))
        
        self.b = self._empty_real((self.mx,self.my))
        self.b1 = self._empty_real((self.mx,self.my))
        self.b2 = self._empty_real((self.mx,self.my))
        self.b3 = self._empty_real((self.mx,self.my))
        self.b4 = self._empty_real((self.mx,self.my))
        
        # for fast transform
        pyfftw.interfaces.cache.enable()
        self.w_to_wh = pyfftw.FFTW(self.w,  self.wh, threads=self.fftw_num_threads,
                                   axes=(-2,-1))
        self.wh_to_w = pyfftw.FFTW(self.wh,  self.w, threads=self.fftw_num_threads,
                                   direction='FFTW_BACKWARD', axes=(-2,-1))
        self.dwhdt_to_dwdt = pyfftw.FFTW(self.dwhdt, self.dwdt, threads=self.fftw_num_threads,
                                         direction='FFTW_BACKWARD', axes=(-2,-1))
        self.u_to_uh = pyfftw.FFTW(self.u,  self.uh, threads=self.fftw_num_threads,
                                   axes=(-2,-1))
        self.uh_to_u = pyfftw.FFTW(self.uh, self.u, threads=self.fftw_num_threads,
                                   direction='FFTW_BACKWARD', axes=(-2,-1))
        self.v_to_vh = pyfftw.FFTW(self.v,  self.vh, threads=self.fftw_num_threads,
                                   axes=(-2,-1))
        self.vh_to_v = pyfftw.FFTW(self.vh, self.v, threads=self.fftw_num_threads,
                                   direction='FFTW_BACKWARD', axes=(-2,-1))
        self.b_to_a = pyfftw.FFTW(self.b, self.a, threads=self.fftw_num_threads,
                                  axes=(-2,-1))
        self.a1_to_b1 = pyfftw.FFTW(self.a1, self.b1, threads=self.fftw_num_threads,
                                    direction='FFTW_BACKWARD', axes=(-2,-1))
        self.a2_to_b2 = pyfftw.FFTW(self.a2, self.b2, threads=self.fftw_num_threads,
                                    direction='FFTW_BACKWARD', axes=(-2,-1))
        self.a3_to_b3 = pyfftw.FFTW(self.a3, self.b3, threads=self.fftw_num_threads,
                                    direction='FFTW_BACKWARD', axes=(-2,-1))
        self.a4_to_b4 = pyfftw.FFTW(self.a4, self.b4, threads=self.fftw_num_threads,
                                    direction='FFTW_BACKWARD', axes=(-2,-1))

        # ṣpectral filter
        try:
            self.fltr
        except AttributeError:
            self._init_filter()
    

    def init_field(self, field="Taylor-Green", t=0.0, kappa=2., delta=0.005, sigma= 15./np.pi):
        """
        Inital the vorticity field. Different fields are hard coded, i.e. Taylor-Green vortex, double shear layer,
        McWilliams random vorticity realisation. User-defined fields can also be passed.

            Params:
                field : string
                    - Type of field to initialise, or actuall field as a numpy array
                t : float
                    - time at which to compute field (onlt for the Taylor-Green solution)
                other : varies
                    - additional parameters, field dependent
        """
        if(type(field)==str):
            if(field=="TG" or field=="Taylor-Green"):
                self.w[:,:] = 2 * kappa * np.cos(kappa * self.x) * np.cos(kappa * self.y[:, np.newaxis]) *\
                              np.exp(-2 * kappa**2 * t / self.Re)
            elif(field=="SL" or field=="Shear Layer"):
                self.w[:,:] = delta * np.cos(self.x) - sigma * np.cosh(sigma * (self.y[:,np.newaxis] -\
                              0.5*np.pi))**(-2)
                self.w [:,:]+= delta * np.cos(self.x) + sigma * np.cosh(sigma * (1.5*np.pi -\
                               self.y[:,np.newaxis]))**(-2)
            elif(field=="McWilliams" or field=="MW84"):
                self.McWilliams1984()
            else:
                print("The specified field type %s is unknown.\nAvailable initial fields are"+\
                      ": \"Taylor-Green\", \"Shear Layer\"." % field)
        elif(type(field)==np.ndarray):
            if(field.shape==(self.nx, self.ny)):
                self.w[:,:] = field
            else:
                print("Specified velocity field does not match grid initialized.")


    # bit-aligned storage arrays for pyFFTW
    def _empty_real(self, *args): 
        shape = (self.nx, self.ny)
        for sp in args:
            shape = sp
        if self.FFTW:
            out = pyfftw.empty_aligned(shape, dtype='float64')
            out.flat[:] = 0.
            return out
        else:
            return np.zeros(shape, dtype='float64')


    def _empty_imag(self, *args):
        shape = (self.nx, self.nk)
        for sp in args:
            shape = sp
        if self.FFTW:
            out = pyfftw.empty_aligned(shape, dtype='complex128')
            out.flat[:] = 0.
            return out
        else:
            return np.zeros(shape, dtype='complex128')


    def get_u(self):
        self.uh[:,:] = self.ky[:,np.newaxis]*self.psih[:, :]
        self.uh_to_u()


    def get_v(self):
        self.vh[:,:] = -self.kx[:self.nk]*self.psih[:, :]
        self.vh_to_v()
        

    def McWilliams1984(self):
        """
        Generates McWilliams vorticity field, see:
            McWilliams (1984), "The emergence of isolated coherent vortices in turbulent flow"
        """
        # generate variable
        self.k2 = self.kx[:self.nk]**2 + self.ky[:,np.newaxis]**2
        self.fk = self.k2 != 0.0

        # ensemble variance proportional to the prescribed scalar wavenumber function
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
        self.w[:,:] = np.fft.irfft2(wh)


    def _init_filter(self):
        """
        Exponential filter, designed to completely dampens higest modes to machine accuracy
        """
        cphi = 0.65*np.max(self.kx)
        wvx = np.sqrt(self.k2)
        filtr = np.exp(-self.filterfac*(wvx-cphi)**4.)
        filtr[wvx<=cphi] = 1.
        self.fltr = filtr


    def _cfl_limit(self):
        """
        Adjust time-step based on the courant condition
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
                - desired order of the method, default is 3rd order
        """
        # iniitalise field
        self.w0[:, :] = self.w[:, :]

        for k in range(s, 0, -1):
            # invert Poisson equation for the stream function (changes to k-space)
            self._get_psih()

            # get convective forces (resets dwhdt)
            self._add_convection()

            # add diffusion (changes to C-space)
            self._add_diffusion()

            # step in time
            self.w[:, :] = self.w0[:, :] + (self.dt/k) * self.dwdt[:, :]

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
            D = 1/Re * (-k_x^2 -k_y^2) * \hat{\omega}
        
        Note: This resets the value in self.w when called
        """
        self.dwhdt[:, :] = self.dwhdt[:, :] - self.ReI*self.k2*self.wh[:, :]
        self.dwhdt_to_dwdt()


    def _add_convection(self):
        """
        Convective term
            N = d/dx \psi * d/dy \omega - d/dy \psi * d/dx \omega
        To prevent alliasing, we zero-pad the array before using the
        convolution theorem to evaluate it in physical space.
        """
        # padded arrays
        j1f_padded = np.zeros((self.mx,self.mk),dtype='complex128')
        j2f_padded = np.zeros((self.mx,self.mk),dtype='complex128')
        j3f_padded = np.zeros((self.mx,self.mk),dtype='complex128')
        j4f_padded = np.zeros((self.mx,self.mk),dtype='complex128')
        
        j1f_padded[self.padder, :self.nk] = 1.0j*self.kx[:self.nk     ]*self.psih[:, :]
        j2f_padded[self.padder, :self.nk] = 1.0j*self.ky[:, np.newaxis]*self.wh[:, :]
        j3f_padded[self.padder, :self.nk] = 1.0j*self.ky[:, np.newaxis]*self.psih[:, :]
        j4f_padded[self.padder, :self.nk] = 1.0j*self.kx[:self.nk     ]*self.wh[:, :]
        
        # ifft
        j1 = self.a1_to_b1(j1f_padded)
        j2 = self.a2_to_b2(j2f_padded)
        j3 = self.a3_to_b3(j3f_padded)
        j4 = self.a4_to_b4(j4f_padded)
        
        jacp = j1*j2 - j3*j4
        
        jacpf = self.b_to_a(jacp)
        
        self.dwhdt[:, :] = jacpf[self.padder, :self.nk]*self.pad**(2) # this term is the result of padding


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


    def tke(self):
        ke = .5*self._spec_variance(np.sqrt(self.k2)*self.psih)
        return ke.sum()


    def enstrophy(self):
        eps = .5*abs(self.w)**2
        return eps.sum(axis=(-2,-1))


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

    
    def plot_spec(self, res=200):
        self._compute_spectrum(200)
        plt.figure(figsize=(6,6))
        plt.plot(self.k, self.E, '-k', label="E(k)")
        plt.xlabel("k")
        plt.ylabel("E(k)")
        plt.legend()
        plt.show()


    def write(self, folder, iter):
        s = np.zeros(self.ny); s[0]=self.time; s[1]=self.dt
        s[2]=self.tke(); s[3]=self.enstrophy()
        np.savetxt(str(folder)+"vort_"+str("%06d"%iter)+".dat", np.vstack((s, self.w)))


    def display(self, complex=False, u_e=None):
        u = self.w
        if complex:
            u = np.real(self.wh)
        if not np.any(u_e)==None:
            u -= u_e
        p=plt.imshow(u, cmap="RdBu")
        plt.colorbar(p)
        plt.xticks([]); plt.yticks([])
        plt.show()


    def display_vel(self):
        if(self.uptodate!=True):
            self.w_to_wh()
            self._get_psih()
            self.get_u()
            self.get_v()
        plt.figure()
        plt.streamplot(self.x, self.y, self.u, self.v)
        plt.xlabel("x"); plt.ylabel("y")
        plt.show()


    def run_live(self, stop, every=100):
        from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
        iterr = 0
        plt.ion()
        fig = plt.figure()
        ax = fig.add_subplot(111)
        im = ax.imshow(self.w, norm=None, cmap="binary_r")
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
                      self.time, stop-self.time, self.tke()))


# if __name__=="__main__":
#     flow = Fluid(128, 128, 1)
#     flow.init_solver()
#     flow.init_field("TG")
#     print(flow.tke())
