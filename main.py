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

if __name__=="__main__":

    # build fluid and solver
    flow = Fluid(512, 512, 2000, pad=1.)
    flow.init_solver()
    flow.init_field("McWilliams")
    
    print("Starting interating on field.\n")
    start_time = t.time()
    finish = 30.0

    # loop to solve
    while(flow.time<=finish):
        flow.update()
        if(flow.it % 3000 == 0):
            print("Iteration \t %d, time \t %f, time remaining \t %f. TKE: %f, ENS: %f" %(flow.it,
                  flow.time, finish-flow.time, flow.tke(), flow.enstrophy()))
            flow.write(folder="Dat/", iter=flow.it/3000)
    # flow.run_live(finish, every=200)

    end_time = t.time()
    print("\nExecution time for %d iterations is %f seconds." %(flow.it, end_time-start_time))
