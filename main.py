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

if __name__=="__main__":

    # build fluid and solver
<<<<<<< HEAD
    flow = Fluid(512, 512, 2000, pad=1.)
    flow.init_solver()
    flow.init_field("McWilliams")
    
    print("Starting interating on field.\n")
    start_time = t.time()
    finish = 30.0
=======
    flow = Fluid(128, 128, 1.e4)
    flow.init_solver()
    q = -np.exp(-((flow.x-np.pi)**2 + (4.0*(flow.y[:,np.newaxis]-np.pi))**2)/(np.pi/3.0)**2)
    # q = np.loadtxt("MW.dat")
    flow.init_field(field=q)

    print("Starting integrating on field.\n")
    start_time = t.time()
    finish = 10.0
>>>>>>> c6ef2771b957183b7abfd84d82b96eb71875def1

    # loop to solve
    while(flow.time<=finish):
        flow.update()
<<<<<<< HEAD
        if(flow.it % 3000 == 0):
            print("Iteration \t %d, time \t %f, time remaining \t %f. TKE: %f, ENS: %f" %(flow.it,
                  flow.time, finish-flow.time, flow.tke(), flow.enstrophy()))
            flow.write(folder="Dat/", iter=flow.it/3000)
    # flow.run_live(finish, every=200)

    end_time = t.time()
    print("\nExecution time for %d iterations is %f seconds." %(flow.it, end_time-start_time))

    flow.display()
=======
        if(flow.it % 500 == 0):
            print("Iteration \t %d, time \t %f, time remaining \t %f. TKE: %f, ENS: %f" %(flow.it,
                    flow.time, finish-flow.time, flow.tke(), flow.enstrophy()))
    #         flow.write(folder="Dat/", iter=flow.it/500)
    # flow.run_live(finish, every=1000)

    end_time = t.time()
    print("\nExecution time for %d iterations is %f seconds." %(flow.it, end_time-start_time))
>>>>>>> c6ef2771b957183b7abfd84d82b96eb71875def1
