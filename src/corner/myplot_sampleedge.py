import numpy as np
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

from matplotlib.animation import FFMpegWriter

def plot_scalar(enum,folder):

    gausspts = 3
    numpts   = 4

    filename = os.path.join(folder, 'build', 'samplevxvert.dat')
    px = np.loadtxt(filename)

    filename = os.path.join(folder, 'build', 'samplevyvert.dat')
    py = np.loadtxt(filename)

    filename = os.path.join(folder, 'build', 'samplevert.dat')
    val = np.loadtxt(filename)

    kk = 14

    ax = plt.axes(projection='3d')
    ax.plot_surface(px[kk].reshape(gausspts, numpts), 
                    py[kk].reshape(gausspts, numpts),
                    val[kk].reshape(gausspts, numpts), cmap='viridis', edgecolor='none')

    plt.show()

    kk = 15

    ax = plt.axes(projection='3d')
    ax.plot_surface(px[kk].reshape(gausspts, numpts), 
                    py[kk].reshape(gausspts, numpts),
                    val[kk].reshape(gausspts, numpts), cmap='viridis', edgecolor='none')

    plt.show()

    for k in range(enum):
        plt.pcolormesh(px[k].reshape(gausspts,numpts), 
                       py[k].reshape(gausspts,numpts), 
                       val[k].reshape(gausspts,numpts), shading='gouraud', cmap='plasma', vmax=2.0, vmin =-0.1)
 
    plt.show() 

    #print(px[k].reshape(gausspts, numpts))
    #print(py[k].reshape(gausspts, numpts))

    #print(val[k].reshape(gausspts, numpts))

