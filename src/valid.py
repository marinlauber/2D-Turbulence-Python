#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Marin Lauber"
__copyright__ = "Copyright 2019, Marin Lauber"
__license__ = "GPL"
__version__ = "1.0.1"
__email__  = "M.Lauber@soton.ac.uk"

import time as t
import numpy as np
from src.fluid import Fluid
from src.field import TaylorGreen,L2,Linf

if __name__=="__main__":

    # build fluid and solver
    flow = Fluid(128, 128, 1.)
    flow.init_solver()
    flow.init_field(TaylorGreen)

    print("Starting integration on field.\n")
    start_time = t.time()
    finish = 0.1

    # loop to solve
    while(flow.time<=finish):

        #  update using RK
        flow.update()

        #  print every 100 iterations
        if (flow.it % 100 == 0):
            print("Iteration \t %d, time \t %f, time remaining \t %f. TKE: %f, ENS: %f" %(flow.it,
                  flow.time, finish-flow.time, flow.tke(), flow.enstrophy()))
    # flow.run_live(finish, every=100)

    end_time = t.time()
    print("\nExecution time for %d iterations is %f seconds." %(flow.it, end_time-start_time))
    
    # get final results
    flow.wh_to_w()
    w_n = flow.w.copy()

    # exact solution
    w_e  = TaylorGreen(flow.x, flow.y,flow.Re, time=flow.time)
    
    # L2-norm and exit 
    print("The L2-norm of the Error in the Taylor-Green vortex on a %dx%d grid is %e." % (flow.nx, flow.ny, L2(w_e - w_n)) )
    print("The Linf-norm of the Error in the Taylor-Green vortex on a %dx%d grid is %e." % (flow.nx, flow.ny, Linf(w_e - w_n)) )
