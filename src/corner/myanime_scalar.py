import numpy as np
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from matplotlib.animation import FFMpegWriter

def anime_scalar(M, N, numsteps, folder, name):

# ---------- Cell grid (x) ----------
    filename = os.path.join(folder, 'build', 'cellgridx.dat')
    cellgx = np.loadtxt(filename)

    # MATLAB: reshape(cellgx, M, N)
    cellgx = cellgx.reshape((M, N))

    # ---------- Cell grid (y) ----------
    filename = os.path.join(folder, 'build', 'cellgridy.dat')
    cellgy = np.loadtxt(filename)

    cellgy = cellgy.reshape((M, N))

    # ============================================================
    # Cell Scalar 
    # ============================================================
    #filename = os.path.join(folder, 'build', f"{name}{mark}.dat")

    #val = np.loadtxt(filename)
    #val = val.reshape(M, N, order='F')
#
#    fig = plt.figure()
#
#    ax = fig.add_subplot(111, projection='3d')
#    surf = ax.plot_surface(cellgx, cellgy, val, cmap='viridis', edgecolor='none')
#
#    fig.colorbar(surf, shrink=0.5)
#    plt.show()

    fig, ax = plt.subplots()
    writer = FFMpegWriter(fps=10, bitrate=1800)

    with writer.saving(fig, f"{name}.mp4", dpi=200):
        for k in range(numsteps):
            kk=k+1
            ax.clear()
            filename = os.path.join(folder, 'build', f"{name}{kk}.dat")
            val = np.loadtxt(filename)
            val = val.reshape(M, N)

            #surf = ax.plot_surface(cellgridx, cellgridy, val, cmap='viridis', edgecolor='none') 
            surf = ax.pcolormesh(cellgx, cellgy, val, cmap='jet', shading='gouraud')
            #cbar = fig.colorbar(surf, ax=ax)
            ax.set_title(f"Step {k}")
            writer.grab_frame()
        fig.colorbar(surf)
