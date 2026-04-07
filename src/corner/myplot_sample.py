import numpy as np
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from matplotlib.animation import FFMpegWriter

def plot_scalar(M, N, folder):

# ---------- Cell grid (x) ----------
    filename = os.path.join(folder, 'build', 'samplevxvert.dat')
    cellgx = np.loadtxt(filename)

    # MATLAB: reshape(cellgx, M, N)
    #cellgx = cellgx.reshape((M, N), order='F')
    cellgx = cellgx.reshape((6*(M+1), N))

    # ---------- Cell grid (y) ----------
    filename = os.path.join(folder, 'build', 'samplevyvert.dat')
    cellgy = np.loadtxt(filename)

    #cellgy = cellgy.reshape((M, N), order='F')
    cellgy = cellgy.reshape((6*(M+1), N))

    print(cellgx)
    #print(cellgy)

    # ============================================================
    # Cell Scalar 
    # ============================================================
#    filename = os.path.join(folder, 'build', f"{name}{mark}.dat")
    filename = os.path.join(folder, 'build', 'samplevert.dat')

    val = np.loadtxt(filename)
    val = val.reshape(6*(M+1), N)

    fig = plt.figure()

    ax = fig.add_subplot(111, projection='3d')
    surf = ax.plot_surface(cellgx, cellgy, val, cmap='viridis', edgecolor='none')

    fig.colorbar(surf, shrink=0.5)
    plt.show()


