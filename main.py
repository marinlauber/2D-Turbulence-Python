#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Marin Lauber"
__copyright__ = "Copyright 2019, Marin Lauber"
__license__ = "GPL"
__version__ = "1.0.1"
__email__  = "M.Lauber@soton.ac.uk"

import time as t
from src.fluid import Fluid
from src.field import DecayingTurbulence

if __name__=="__main__":

    # build fluid and solver
    flow = Fluid(512, 512, 2000, pad=1.)
    flow.init_solver()
    flow.init_field(DecayingTurbulence)

    print("Starting integration on field.\n")
    start_time = t.time()
    finish = 0.01

    # loop to solve
    while(flow.time<=finish):
        flow.update()
        if(flow.it % 100 == 0):
            print("Iteration \t %d, time \t %f, time remaining \t %f. TKE: %f, ENS: %f" %(flow.it,
                  flow.time, finish-flow.time, flow.tke(), flow.enstrophy()))
            flow.write(file="fluid")
    
    # flow.run_live(finish, every=200)

    end_time = t.time()
    print("\nExecution time for %d iterations is %f seconds." %(flow.it, end_time-start_time))