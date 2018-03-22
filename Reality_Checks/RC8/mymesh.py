import numpy as np
from mpl_toolkits.mplot3d import axes3d, Axes3D
import matplotlib.pyplot as plt


# 3-D plot of 2D array w
def mymesh(xvals, yvals, w, xlabel='', ylabel='', zlabel=''):
    x, y = np.meshgrid(xvals, yvals)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_wireframe(x, y, w, rstride=1)
    ax.set_xlabel(xlabel);
    ax.set_ylabel(ylabel);
    ax.set_zlabel(zlabel)
    plt.show()
