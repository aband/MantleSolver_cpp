import numpy as np
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from matplotlib.animation import FFMpegWriter

def plot_single(cellnum, m, n, folder):
    filename = os.path.join(folder, 'build', 'patchx.dat')
    px       = np.loadtxt(filename)

    filename = os.path.join(folder, 'build', 'patchy.dat')
    py       = np.loadtxt(filename)

    filename = os.path.join(folder, 'build', 'patch.dat')
    val      = np.loadtxt(filename)

#    newpx = np.arrange(cellnum*m*n)
#    newpy = np.arrange(cellnum*m*n)

#    for jj in range(cellnum):
#        for nn in range(n):
#            for mm in range(m):
#                newpx[jjkkkkkk] = 

#    fig=plt.figure()

    ax = plt.axes(projection='3d')
    surf = ax.plot_surface(px[11].reshape(m,n), 
                           py[11].reshape(m,n),
                           val[11].reshape(m,n), cmap='viridis', edgecolor='none')

    plt.show() 

    ax = plt.axes(projection='3d')
    surf = ax.plot_surface(px[12].reshape(m,n), 
                           py[12].reshape(m,n),
                           val[12].reshape(m,n), cmap='viridis', edgecolor='none')

    plt.show() 

    ax = plt.axes(projection='3d')
    surf = ax.plot_surface(px[13].reshape(m,n), 
                           py[13].reshape(m,n),
                           val[13].reshape(m,n), cmap='viridis', edgecolor='none')

    plt.show() 

    for k in range(cellnum):
        #ax = fig.add_subplot(111, projection='3d')
        #surf = ax.plot_surface(px[k].reshape(m,n), 
        #                       py[k].reshape(m,n),
        #                       val[k].reshape(m,n), cmap='viridis', edgecolor='none')
        plt.pcolormesh(px[k].reshape(m,n), py[k].reshape(m,n), val[k].reshape(m,n), shading='gouraud', cmap='plasma', vmax=2.0, vmin=-0.1)
    #plt.colorbar(label='Intensity')
    #plt.colorbar()
    plt.show()
