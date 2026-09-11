import numpy as np
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from matplotlib.animation import FFMpegWriter


def plot_edge_scalar(M, N, mark, folder, name):

    # ---------- Vertical Gauss grid (x) ----------
    filename = os.path.join(folder, 'build', 'vertgaussgridx.dat')
    vertgx = np.loadtxt(filename)

    # MATLAB: reshape(vertgx, 3, N*(M+1))
    vertgx = vertgx.reshape((3, N * (M + 1)), order='F')

    revertgx = vertgx[:, :M + 1]
    for i in range(1, N):
        block = vertgx[:, i * (M + 1):(i + 1) * (M + 1)]
        revertgx = np.vstack((revertgx, block))

    # ---------- Vertical Gauss grid (y) ----------
    filename = os.path.join(folder, 'build', 'vertgaussgridy.dat')
    vertgy = np.loadtxt(filename)

    vertgy = vertgy.reshape((3, N * (M + 1)), order='F')

    revertgy = vertgy[:, :M + 1]
    for i in range(1, N):
        block = vertgy[:, i * (M + 1):(i + 1) * (M + 1)]
        revertgy = np.vstack((revertgy, block))

    # ---------- Horizontal Gauss grid (x) ----------
    filename = os.path.join(folder, 'build', 'horigaussgridx.dat')
    horigx = np.loadtxt(filename)

    # MATLAB: reshape(horigx, M*3, N+1)
    horigx = horigx.reshape((M * 3, N + 1), order='F')

    # ---------- Horizontal Gauss grid (y) ----------
    filename = os.path.join(folder, 'build', 'horigaussgridy.dat')
    horigy = np.loadtxt(filename)

    horigy = horigy.reshape((M * 3, N + 1), order='F')

    # ---------- Degrees of freedom ----------
    vertdof = N * (M + 1) * 3
    horidof = M * (N + 1) * 3

    filename = os.path.join(folder, 'build', f"{name}{mark}.dat")
    val = np.loadtxt(filename)

    # Vertical values
    valvert = val[:vertdof]
    valvert = valvert.reshape((3, N * (M + 1)), order='F')

    revalverty = valvert[:, :M + 1]
    for i in range(1, N):
        block = valvert[:, i * (M + 1):(i + 1) * (M + 1)]
        revalverty = np.vstack((revalverty, block))

    # Horizontal values
    valhori = val[vertdof:vertdof + horidof]
    valhori = valhori.reshape((M * 3, N + 1), order='F')

    # plot
    fig, ax = plt.subplots()

    #ax = fig.add_subplot(111, projection='3d')
    #surf = ax.pcolormesh(vertgx, vertgy, valvert, cmap='viridis', shading='gouraud')
    surf = ax.pcolormesh(horigx, horigy, valhori, cmap='viridis', shading='gouraud', vmin=0.0, vmax=0.4)

    fig.colorbar(surf, shrink=0.5)
    plt.show()


