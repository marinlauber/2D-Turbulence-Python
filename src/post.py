#!/usr/bin/env python3
# -*- coding: utf-8 -*-

__author__ = "Marin Lauber"
__copyright__ = "Copyright 2019, Marin Lauber"
__license__ = "GPL"
__version__ = "1.0.1"
__email__  = "M.Lauber@soton.ac.uk"

import numpy as np
import matplotlib.pyplot as plt
import os

if __name__=="__main__":

    print("Cleaning workspace..")
    os.system("rm *.png vid.mp4")
    print("Done.")

    print("Generating .png files and assembling movie...")
    i = 0

    while(os.path.isfile("vort_"+"%06d"%(i+1)+".dat")):

        w = np.genfromtxt("vort_"+"%06d"%(i+1)+".dat", delimiter=" ")
        time=w[0,0]; dt=[0,1]; w=np.round(w[1:,:],10)
        if i==0: min_=np.min(w);max_=np.max(w)

        plt.figure(figsize=(8,6))
        plt.contourf(np.linspace(0, 1, w.shape[0]), np.linspace(0, 1, w.shape[1]),
                     np.round(w, 5), cmap="RdBu",levels=51, vmin=min_, vmax=max_)
        plt.colorbar()
        plt.xticks([]); plt.yticks([])
        plt.title(f"Time: {time:.2f} s")
        plt.savefig("vort_"+"%06d"%(i+1)+".png", dpi=600);i+=1
        plt.close()

    os.system("ffmpeg -r 24 -i vort_%06d.png -vcodec libx264 -crf 25 -pix_fmt yuv420p vid.mp4")
    os.system("rm *.png")
    print("Done.\nExiting.")
