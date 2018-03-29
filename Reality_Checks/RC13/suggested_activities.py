import numpy as np
from get_potential_energy import *
import scipy.optimize as opt
# from nelder_mead import *


def sa1():
    # ig = initial guess, rad = radius, k = # of iterations
    z2 = 1
    origin = np.array([0, 0, 0])
    v1 = np.array([0, 0, z2])
    # ig = np.array([origin, v1, [1, 0, 0], [0, 1, 0], [1, 1, 0]])
    ig = np.array([1., 2., 3., 4., 5.])
    rad = 1
    k = 50
    xmin = opt.minimize(get_potential_energy, ig, method='Nelder-Mead', tol=.0005)
    # xmin = neldermead(get_potential_energy, ig, rad, k)
    # print('A minimum of %.6f' % get_potential_energy(np.mean(xmin, axis=1)), 'occurs at (%.6f, %.6f)'
    #       % (np.mean(xmin, axis=1)[0], np.mean(xmin, axis=1)[1]))
    print(xmin)
    return


def sa2():
    return


def sa3():
    return


def sa4():
    return


def sa5():
    return


def sa6():
    return


def sa7():
    return
