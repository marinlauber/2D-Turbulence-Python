#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Marin Lauber"
__copyright__ = "Copyright 2019, Marin Lauber"
__license__ = "GPL"
__version__ = "1.0.1"
__email__  = "M.Lauber@soton.ac.uk"

"""
Change log:
    16.10.2019
    do: cfl_limit() needs to be corrected
"""

import numpy as np
import matplotlib.pyplot as plt

class fluid(object):

    def __init__(self, nx, ny, Re, dt=0.0001, pad=3./2.):

        # input data
        self.nx = nx
        self.ny = ny; self.nk = self.ny//2+1
        self.Re = Re; self.ReI = 1./self.Re
        self.dt = dt
        self.pad = pad
        self.time = 0.
        self.uptodate = False

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

    
    def _init_field(self, field="Taylor-Green", t=0.0, kappa=2., delta=0.005, sigma= 15./np.pi):
        """
        Inital flow field following the Taylor-Green solution of the Navier-Stokes or a double shear layer.

            Params:
                field: string
                    -Type of field to initialise
        """
        if(field=="TG" or field=="Taylor-Green"):
            self.w = 2 * kappa * np.cos(kappa * self.x) * np.cos(kappa * self.y[:, np.newaxis]) *\
                     np.exp(-2 * kappa**2 * t / self.Re)
        elif(field=="SL" or field=="Shear Layer"):
            self.w = delta * np.cos(self.x) - sigma * np.cosh(sigma * (self.y[:,np.newaxis] - 0.5*np.pi))**(-2)
            self.w += delta * np.cos(self.x) + sigma * np.cosh(sigma * (1.5*np.pi - self.y[:,np.newaxis]))**(-2)
        else:
            print("The specified field type %s is unknown."&
                +"\nAvailable initial fields are: \"Taylor-Green\", \"Shear Layer\"." % field)


    def _init_solver(self):

        # initialise array required for solving
        self.wh = np.empty((self.nx,self.nk), dtype=np.complex128)
        self.w0 = np.empty((self.nx,self.nk), dtype=np.complex128)
        self.psih = np.zeros((self.nx,self.nk), dtype=np.complex128)
        self.dwhdt = np.zeros((self.nx,self.nk), dtype=np.complex128)

        # utils
        self.k2 = self.kx[:self.nk]**2 + self.ky[:,np.newaxis]**2
        self.fk = self.k2 != 0.0
        self.mx = int(self.pad * self.nx)
        self.my = int(self.pad * self.nk)
        self.padder = np.ones(self.mx, dtype=bool)
        self.padder[int(self.nx/2):int(self.nx*(self.pad-0.5)):] = False

        # populate those arrays
        self.w_to_wh()
        self.get_psih()
        self.w0 = self.wh


    def w_to_wh(self):
        self.wh = np.fft.rfft2(self.w, axes=(-2,-1))
        

    def wh_to_w(self):
        self.w = np.fft.irfft2(self.wh, axes=(-2,-1))


    def get_u(self):
        """
        Spectral differentitation to get:
            u = d/dy \psi
        """
        self.u = np.fft.irfft2(self.ky[:,np.newaxis]*self.psih)


    def get_v(self):
        """
        Spectral differentitation to get:
            v = -d/dx \psi
        """
        self.v = -np.fft.irfft2(self.kx[:self.nk]*self.psih)


    def cfl_limit(self):
        """
        Adjust time-step based on the courant condition
        
        Note: this assumes that you initial velocity field is correctly normalized.
        """
        self.get_u()
        self.get_v()
        self.dt = 0.001*np.min(np.hstack((self.dx, self.dy)))/np.max(np.hstack((self.u, self.v)))


    def get_psih(self):
        """
        Spectral stream-function from spectral vorticity
            hat{\psi} = \hat{\omega} / (k_x^2 + k_y^2)
        """
        self.psih[self.fk] = self.wh[self.fk] / self.k2[self.fk]
    

    def update(self, s=3):
        """
        Low-storage S-order Runge-Kutta method from Jameson, Schmidt and Turkel (1981)
        Input:
            s : float
                desired order of the method, default is 3rd order
        """
        # iniitalise field
        self.w0 = self.wh

        for k in range(s,0,-1):

            # invert Poisson equation for the stream function
            self.get_psih()

            # get convective forces (reset dwhdt)
            self.add_convection()

            # add diffusion
            self.add_diffusion()

            # step in time
            self.wh = self.w0 + (self.dt/k) * self.dwhdt

        self.time += self.dt
        self.cfl_limit()


    def add_diffusion(self):
        """
        Diffusion term of the Navier-Stokes
            D = p * 1/Re * (-k_x^2 -k_y^2) * \hat{\omega}
        
        Note: This resets the value in self.w when called
              The penalty value is required for the RK method
        """
        self.dwhdt -= self.ReI*self.k2*self.wh


    def add_convection(self):
        """
        Convective term
            N = -d/dy \psi * d/dx \omega + d/dx \psi * d/dy \omega
        To prevent alliasing, we zero-pad the array before using the
        convolution theorem to evaluate it in physical space.
        """
        self.dwhdt = -1*self.convolve(1j*self.ky[:,np.newaxis]*self.psih, 1j*self.kx[:self.nk]*self.wh)
        self.dwhdt += self.convolve(1j*self.kx[:self.nk]*self.psih, 1j*self.ky[:,np.newaxis]*self.wh)


    def convolve(self, a, b):
        """
        Evaluate convolution sum. This involves three transforms
        """
        # zero-padded temp arrays
        tmp = np.zeros((self.mx,self.my), dtype=np.complex128)

        # padd input arrays
        tmp[self.padder,:self.nk] = a
    
        # fft with these new coeff, padded with zeros
        r = np.fft.irfft2(tmp, axes=(-2,-1))

        # padd input arrays
        tmp[self.padder,:self.nk] = b

        # multiplication in physical space, this saves one temp array
        r *= np.fft.irfft2(tmp, axes=(-2,-1))*self.pad**(2)

        # transform back to k-space, need normalisation, re-assigne to old array
        tmp = np.fft.rfft2(r, axes=(-2,-1))
        
        return tmp[self.padder, :self.nk] # truncate fourier modes


    def write_vort(self):
        np.savetxt("vort_"+str(self.dt)+".dat", self.w)


    def show_vort(self):
        plt.figure()
        plt.imshow(self.w)
        plt.show()
    

    def show_vel(self):
        if(self.uptodate!=True):
            self.w_to_wh()
            self.get_psih()
            self.get_u()
            self.get_v()
        plt.figure()
        plt.quiver(self.x, self.y, self.u, self.v)
        plt.xlabel("x"); plt.ylabel("y")
        plt.show()
        

# if __name__=="__main__":
#     flow = fluid(128, 128, 1)
