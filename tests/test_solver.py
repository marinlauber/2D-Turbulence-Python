#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Marin Lauber"
__copyright__ = "Copyright 2019, Marin Lauber"
__license__ = "GPL"
__version__ = "1.0.1"
__email__  = "M.Lauber@soton.ac.uk"

import numpy as np
import matplotlib.pyplot as plt
from src.fluid import Fluid

def test_diagnostics():
    flow = Fluid(64, 64, 1.)
    flow.init_solver()
    flow.init_field("McWilliams")
    flow._get_psih()
    assert np.isclose(np.round(flow.tke(),1), 0.5, atol=1e-3),\
           "Error: TKE do not match"
    assert flow.enstrophy!=0.0, "Error: Enstrophy is zero"
    flow._compute_spectrum(200)


def test_plots():
    flow = Fluid(16, 16, 1.)
    flow.init_solver()
    flow.init_field("Shear Layer")
    plt.ion()
    flow.plot_spec()
    plt.close()
    flow.display()
    plt.close()
    flow.display_vel()
    plt.close()


def test_update():
    # build field
    flow = Fluid(32, 32, 1.)
    flow.init_solver()
    flow.init_field("Taylor-Green")
    # start update
    while(flow.time<=0.1):
        flow.update()
    # get final results
    w_n = flow.w.copy()
    # exact solution
    flow.init_field("Taylor-Green", t=flow.time)
    assert np.allclose(w_n, flow.w, atol=1e-6), "Error: solver diverged."


# if __name__=="__main__":
#     test_diagnostics()
#     test_plots()
#     test_update()
