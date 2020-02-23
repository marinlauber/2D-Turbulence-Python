#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Marin Lauber"
__copyright__ = "Copyright 2019, Marin Lauber"
__license__ = "GPL"
__version__ = "1.0.1"
__email__  = "M.Lauber@soton.ac.uk"

import time as t
import numpy as np
from fluid import Fluid
import matplotlib.pyplot as plt

def f(x, y, d=0.1, coo=True): 
    x_0 = abs((x[0]+x[-1])*.5)
    y_0 = abs((y[0]+y[-1])*.5)
    x = x-x_0
    y = y-y_0
    x = x/x.max()/2.5
    y = y/y.max()/2.5
    q =  np.exp(-((x-d)**2 + (y[:, np.newaxis])**2)/(.1*d)) +\
               (2*int(coo)-1)*np.exp(-((x+d)**2 + (y[:, np.newaxis])**2)/(.1*d))
    return q

def jet(x, y, delta=0.5*np.pi):
    u = np.where(abs(y[:,np.newaxis]-np.pi)<delta, 1, 0) + 0.*x
    v = np.zeros_like(u)
    return u, v


def psi(x, y, R=0.16*np.pi, eps=0.1):
    return eps/(2*np.pi)*np.sqrt(np.exp(1-(((x-np.pi)**2 + (y[:,np.newaxis]-np.pi)**2)/R**2)))

def vortex_conv(x, y, uinf=1, R=0.16*np.pi, eps=0.1):
    psi = eps/(2*np.pi)*np.sqrt(np.exp(1-(((x-np.pi)**2 + (y[:,np.newaxis]-np.pi)**2)/R**2)))
    u = uinf + (y[:,np.newaxis]/R)*psi
    v = -(x/R)*psi
    return u, v

if __name__=="__main__":

    # build fluid and solver
    flow = Fluid(64, 64, 1.e2)
    flow.init_solver()
    # q = -np.exp(-((flow.x-np.pi)**2 + (4.0*(flow.y[:,np.newaxis]-np.pi))**2)/(np.pi/3.0)**2)
    # q = np.loadtxt("MW.dat")
    u, v, = vortex_conv(flow.x,flow.y)
    q = flow.vorticity(u, v)
    flow.init_field(field=q)
    # flow.display()

    print("Starting integrating on field.\n")
    start_time = t.time()
    finish = 5.0

    # loop to solve
    # while(flow.time<=finish):
    #     flow.update()
    #     if(flow.it % 500 == 0):
    #         print("Iteration \t %d, time \t %f, time remaining \t %f. TKE: %f, ENS: %f" %(flow.it,
    #                 flow.time, finish-flow.time, flow.tke(), flow.enstrophy()))
    # #         flow.write(folder="Dat/", iter=flow.it/500)
    flow.run_live(finish, every=50)

    end_time = t.time()
    print("\nExecution time for %d iterations is %f seconds." %(flow.it, end_time-start_time))
