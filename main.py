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

def merger(x, y, d=0.1, coo=True): 
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
    flow = Fluid(256, 256, 1e12)
    flow.init_field(field=merger(flow.x, flow.y))
    flow.init_solver()

    print("Starting interating on field.\n")
    start_time = t.time()
    iterr = 0
    finish = 20.0

    # loop to solve
    while(flow.time<=finish):
        flow.update()
        iterr += 1
        if(iterr % 500 == 0):
            print("Iteration \t %d, time \t %f, time remaining \t %f" %(iterr,
                                                                        flow.time,
                                                                        finish-flow.time))
            flow.write("Dat/", iterr/500)
    # flow.run_live(finish, every=1000)

    end_time = t.time()
    print("\nExecution time for %d iterations is %f seconds." %(iterr, end_time-start_time))

    flow.show_vort()
