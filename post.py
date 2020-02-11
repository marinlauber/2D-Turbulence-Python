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

def plot(x, y, ID):
    plt.plot(x, y, '-k', label=ID)
    plt.legend()
    plt.savefig(ID+".png")
    plt.close()


def save_image(data, fn):
    sizes = np.shape(data)
    height = float(sizes[0])
    width = float(sizes[1])
    
    fig = plt.figure()
    fig.set_size_inches(width/height, 1, forward=False)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    ax.imshow(data, cmap="RdBu")
    plt.savefig(fn, dpi = height) 
    plt.close()


def save_contour(data, fn):
    plt.figure(figsize=(8,6))
    plt.contourf(np.linspace(0, 1, data.shape[0]), np.linspace(0, 1, data.shape[1]),
                    np.round(data, 5), cmap="RdBu",levels=51)
    plt.colorbar()
    plt.xticks([]); plt.yticks([])
    plt.title(f"Time: {time:.2f} s")
    plt.savefig(fn, dpi=600)
    plt.close()


if __name__=="__main__":

    print("Cleaning workspace..")
    os.system("rm *.png vid.mp4")
    print("Done.")

    print("Generating .png files and assembling movie...")
    i = 0

    t=[]; dt=[]; tke=[]; ens=[]

    while(os.path.isfile("vort_"+"%06d"%(i+1)+".dat")):

        w = np.genfromtxt("vort_"+"%06d"%(i+1)+".dat", delimiter=" ")
        time = w[0,0]; t.append(time); 
        dt.append(w[0,1]); tke.append(w[0,2]); ens.append(w[0,3])
        w=np.round(w[1:,:],10)

        save_image(w, "vort_"+"%06d"%(i+1)+".png")
        i += 1

    os.system("ffmpeg -r 24 -i vort_%06d.png -vcodec libx264 -crf 25 -pix_fmt yuv420p vid.mp4")
    os.system("rm *.png")

    plot(t, dt, "dt")
    plot(t, tke, "Tke")
    plot(t, ens, "Enstrophy")

    print("Done.\nExiting.")
