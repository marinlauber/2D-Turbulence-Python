#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Marin Lauber"
__copyright__ = "Copyright 2020, Marin Lauber"
__license__ = "GPL"
__version__ = "1.0.1"
__email__  = "M.Lauber@soton.ac.uk"

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import os
from math import ceil as ceil
import netCDF4 as nc

try:
    plt.style.use('mystyle')
except OSError:
    print("Using default ploting scheme...")

def plot(x, y, ID):
    plt.plot(x, y, '-k', label=ID)
    plt.legend()
    plt.savefig(ID+".png")
    plt.close()


def save_image(data, fn, cm="RdBu"):
    sizes = np.shape(data)
    height = float(sizes[0]); width = float(sizes[1])
    fig = plt.figure()
    fig.set_size_inches(width/height, 1, forward=False)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    ax.contourf(data, cmap=cm)
    plt.savefig(fn, dpi = 4*height, bbox_inches="tight", pad_inches=0) 
    plt.close()


def save_contour(data, fn, time, cm="RdBu"):
    plt.figure(figsize=(8,6))
    plt.contourf(np.linspace(0, 1, data.shape[0]), np.linspace(0, 1, data.shape[1]),
                 np.round(data, 5), cmap=cm,levels=51)
    plt.colorbar()
    plt.xticks([]); plt.yticks([])
    plt.title(f"Time: {time:.2f} s")
    plt.savefig(fn, dpi=600)
    plt.close()


def save_comp(w, fn, time, cm="PRGn", res=(1920,1080)):
    nx = w.shape[0]; nk=w.shape[1]//2+1
    kx, kk = _wavenumber(nx, nk)
    wh = np.fft.rfft2(w, axes=(-2,-1))
    k2 = kx[:,np.newaxis]**2+kk**2
    psih = get_psi(wh, k2)
    # u, v = get_velocity(psih, kx, kk)
    kE, k, E = get_tke(psih, k2)
    k, O = get_ens(wh, k2)

    fig = plt.figure(figsize=(12,6.75))
    ax1 = plt.subplot2grid((2,2), (0,0), rowspan=2)
    ax2 = plt.subplot2grid((2,2), (0,1))
    ax3 = plt.subplot2grid((2,2), (1,1))

    # plot vorticity on two left cells
    p=ax1.imshow(np.round(w, 5),cmap=cm,extent=[0,1,0,1])
    divider = make_axes_locatable(ax1)
    cax = divider.append_axes('right', size='5%', pad=0.15)
    fig.colorbar(p, cax=cax, orientation='vertical')

    # set labels
    ax1.set_xlabel(r'$x/2\pi$'); ax1.set_ylabel(r'$y/2\pi$')
    ax1.set_title(r'$\omega:=\nabla^2\psi$')

    # plot TKE on top right
    ax2.loglog(k, k*E/np.sum(E), '-k')
    ax2.loglog([8,80],[0.08,0.000008], ':k', alpha=0.5, lw=0.5)
    ax2.text(92,0.000004,r'$k^{-4}$')
    ax2.set_xlabel('Wavenumber'); ax2.set_ylabel(r'$E(k)/\sum E(k)$')
    ax2.set_title('Energy Spectrum')
    ax2.set_ylim(1e-8,1)

    # plot ENS on bottom right
    ax3.loglog(k, k*O/np.sum(O), '-k')
    ax3.loglog([8,80],[0.08,0.8*10**(-5/3.)], ':k', alpha=0.5, lw=0.5)
    ax3.text(92,0.6*10**(-5/3.),r'$k^{-5/3}$')
    ax3.set_ylim(1e-5,1)
    ax3.set_xlabel('Wavenumber'); ax3.set_ylabel(r'$S(k)/\sum S(k)$')
    ax3.set_title('Enstrophy Spectrum')

    # finalize
    fig.suptitle(f"Time: {time:.2f} s",x=0.25,y=0.98,fontsize=16)
    plt.tight_layout()
    plt.savefig(fn, dpi=ceil(res[0]/12.))
    plt.close()


def _wavenumber(*args):
    k_i = ()
    for arg in args:
        k = np.fft.fftfreq(arg, d=1./arg)
        k_i = k_i + (k,)
    return k_i


def get_psi(wh, k2):
    k2I = np.zeros_like(k2)
    fk = k2 != 0.0
    k2I[fk] = 1./k2[fk]
    psih = wh * k2I
    return psih


def get_velocity(psih, kx, ky):
    uh =  1j*ky[:,np.newaxis]*psih
    vh = -1j*kx*psih
    return np.fft.irfft2(uh, axes=(-2,-1)), np.fft.irfft2(vh, axes=(-2,-1))


def get_tke(psih, k2, res=256):
    tke = np.real(0.5*k2*psih*np.conj(psih))
    kmod = np.sqrt(k2)
    k = np.arange(1, k2.shape[0], 1, dtype=np.float64) # nyquist limit for this grid
    E = np.zeros_like(k)
    dk = (np.max(k)-np.min(k))/res
    for i in range(len(k)):
        E[i] += np.sum(tke[(kmod<k[i]+dk) & (kmod>=k[i]-dk)])
    kE = np.sum(k*E)/np.sum(E)
    return kE, k, E


def get_ens(wh, k2, res=256):
    ens = np.real(0.5*k2*wh*np.conj(wh))
    kmod = np.sqrt(k2)
    k = np.arange(1, k2.shape[0], 1, dtype=np.float64) # nyquist limit for this grid
    O = np.zeros_like(k)
    dk = (np.max(k)-np.min(k))/res
    for i in range(len(k)):
        O[i] += np.sum(ens[(kmod<k[i]+dk) & (kmod>=k[i]-dk)])
    return k, O


if __name__=="__main__":
    
    print("Cleaning workspace..")
    os.system("rm -f *.png vid.mp4")
    print("Done.")

    print("Generating .png files and assembling movie...")

    # the resolution we want
    res=(2560,1440)

    # load the data
    data = nc.Dataset("fluid.nc")

    # save avery snapshot
    for i in range(len(data["t"])):
        w = data["w"][i,:,:]
        save_comp(w, "vort_"+"%06d"%(i+1)+".png", time=data["t"][i], cm="RdBu", res=res)
        i += 1

    # generate the movie
    os.system("ffmpeg -r 24 -i vort_%06d.png -s "+str(res[0])+"x"+str(res[1])+" -vcodec libx264 -crf 25 -pix_fmt yuv420p vid.mp4")
    os.system("rm *.png")

    print("Done.\nExiting.")
