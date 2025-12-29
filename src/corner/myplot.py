import numpy as np
import os
import matplotlib.pyplot as plt

def plot_velocity(M, N, mark, folder, name):
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

    # ---------- Velocity x-component filename ----------
    filename = os.path.join(
        folder, 'build', f"{name}x{mark}.dat"
    )

    # ============================================================
    # X component
    # ============================================================
    filename = os.path.join(folder, 'build', f"{name}x{mark}.dat")
    val = np.loadtxt(filename)

    # Vertical values
    valvert = val[:vertdof]
    valvert = valvert.reshape((3, N * (M + 1)), order='F')

    revalvertx = valvert[:, :M + 1]
    for i in range(1, N):
        block = valvert[:, i * (M + 1):(i + 1) * (M + 1)]
        revalvertx = np.vstack((revalvertx, block))

    # Horizontal values
    valhorix = val[vertdof:vertdof + horidof]
    valhorix = valhorix.reshape((M * 3, N + 1), order='F')

    # ============================================================
    # Y component
    # ============================================================
    filename = os.path.join(folder, 'build', f"{name}y{mark}.dat")
    val = np.loadtxt(filename)

    # Vertical values
    valvert = val[:vertdof]
    valvert = valvert.reshape((3, N * (M + 1)), order='F')

    revalverty = valvert[:, :M + 1]
    for i in range(1, N):
        block = valvert[:, i * (M + 1):(i + 1) * (M + 1)]
        revalverty = np.vstack((revalverty, block))

    # Horizontal values
    valhoriy = val[vertdof:vertdof + horidof]
    valhoriy = valhoriy.reshape((M * 3, N + 1), order='F')

    # ============================================================
    # Plot: Vertical Gauss Points
    # ============================================================
#    plt.figure()
#    plt.quiver(revertgx, revertgy, revalvertx, revalverty)
#    plt.streamplot(
#        revertgx, revertgy,
#        revalvertx, revalverty,
#        linewidth=2, color='k', density=1
#    )
#    plt.title('Vertical Gauss Points')

    # ============================================================
    # Plot: Horizontal Gauss Points
    # ============================================================
#    plt.figure()
#    plt.quiver(horigx, horigy, valhorix, valhoriy)
#    plt.streamplot(
#        horigx.T, horigy.T,
#        valhorix.T, valhoriy.T,
#        linewidth=2, color='k', density=1
#    )
#    plt.title('Horizontal Gauss Points')

    # ============================================================
    # Plot: All Gauss Points Together
    # ============================================================
    plt.figure()
    plt.quiver(revertgx, revertgy, revalvertx, revalverty)
    plt.quiver(horigx, horigy, valhorix, valhoriy)
    plt.title('Plot Velocities on Gauss Points All together')

    plt.show()

