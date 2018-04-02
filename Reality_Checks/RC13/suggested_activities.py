import numpy as np
from get_potential_energy import *
import scipy.optimize as opt
from nelder_mead import *


def sa1():
    """
    Write a function that returns the potential energy U=∑i<j1ri j12-1ri j1
    where ri_j is given at the top of p. 581. Apply Nelder–Mead to find the
    minimum energy for n=5. Try several initial guesses until you are convinced
    you have the absolute minimum. How many steps are required?To help you
    check your potential function when n=5, here is one correct input, output
    pair. U(1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0) = -6.0102023319615911
    """
    # fixed points in our model
    origin = np.array([0, 0, 0])
    z2 = 1
    v1 = np.array([0, 0, z2])

    # ig = np.array([origin, v1, [1, 0, 0], [0, 1, 0], [1, 1, 0]])
    ig = np.array([origin, v1])

    # need to first convert the multidimensional array into a vector array
    m, n = np.shape(ig)
    x = m * n
    ig = np.reshape(ig, x)

    # Chilton advised we should use scipy.optimize over his code
    x_min = opt.minimize(get_potential_energy, ig, method='Nelder-Mead')

    print(x_min)
    return


def sa2():
    """
    Plot the five atoms in the minimum energy configuration as circles,
    and connect all circles with line segments to view the conformed
    molecule. The line segments should form triangles and there should
    be three or four lines connecting each atom. You are welcome to use
    Python or Mathematica.
    """
    return


def sa3():
    """
    Write a function that returns the gradient of U. Apply a Python
    minimization function that uses the function and the gradient for
    the n=5 case. Find the minimum energy as before.
    To help you check your gradient function when n=5, here is one
    correct input, output pair.
    ∇U(1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0) = [0.65625, 0.0,
    0.65625, 0.65625,0.65625,-1.3125, 0.79891,-1.45516, 0.79891,-1.3125,
    0.65625, 0.65625,-0.79891, 0.14266,-0.79891]
    """
    return


def sa4():
    """
    Use one of the functions in SciPy Optimization to find the global minimum
    of your potential function when n=5 using only the potential function itself(
    not the gradient). You cannot use Nelder-Mead for this.
    """
    return


def sa5():
    """
    Apply the methods used in (1), (3), and (4) when n=6. Rank the methods
    according to reliability and efficiency.
    """
    return


def sa6():
    """
    Plot the six atoms in the minimum energy configuration as circles, and
    connect all circles with line segments to view the conformed molecule.
    The line segments should form triangles and there should be four lines
    connecting each atom. You are welcome to use Python or Mathematica.
    """
    return


def sa7():
    """
    Determine and plot minimum-energy conformations for larger n. Information
    on minimum-energy Lennard-Jones clusters for n up to several hundred is posted
    at the link provided in (b) above, so your answers can be readily checked. You
    should do at least n=8.
    """
    return
