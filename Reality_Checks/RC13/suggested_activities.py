import numpy as np
from get_potential_energy import *
import scipy.optimize as opt
from nelder_mead import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def sa1():
    """
    Write a function that returns the potential energy U=∑i<j (1/r_ij^12 -1/r_ij^6)
    where r_ij is given at the top of p. 581. Apply Nelder–Mead to find the
    minimum energy for n=5. Try several initial guesses until you are convinced
    you have the absolute minimum. How many steps are required? To help you
    check your potential function when n=5, here is one correct input, output
    pair. U(1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0) = -6.0102023319615911
    """
    # fixed points in our model
    origin = np.array([0, 0, 0])
    z2 = 1
    v1 = np.array([0, 0, z2])

    ig = np.array([origin, v1, [1, 0, 0], [0, 1, 0], [1, 1, 0]])
    # ig = np.array([origin, v1])

    # need to first convert the multidimensional array into a vector array
    m, n = np.shape(ig)
    x = m * n
    ig = np.reshape(ig, x)

    # Chilton advised we should use scipy.optimize over his code
    x_min = opt.minimize(U, ig, method='Nelder-Mead')

    # check
    print("Value should be -6.0102023319615911: ",
          U([1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0]))

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
    # optimal node location for n = 5
    # http://doye.chem.ox.ac.uk/jon/structures/LJ/tables.150.html
    x = np.array([-0.26047, 0.26047, -0.41449, -0.19441, 0.60890])
    y = np.array([0.73631, -0.73632, -0.36526, 0.28435, 0.08091])
    z = np.array([0.47271, -0.47271, 0.34056, -0.55004, 0.20949])
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # plot each point
    ax.scatter(x, y, z, 'bo')

    # plot line connecting the points
    # https://stackoverflow.com/questions/41563347/connecting-points-to-a-central-point-on-3d-scatter-python?rq=1
    for a, b, c in zip(x, y, z):
        for d, e, f in zip(x, y, z):
            ax.plot3D([a, d], [b, e], [c, f], 'b-')

    plt.title('Approximate Solution for $n=5$', fontsize=24)
    plt.savefig('sa2_fig.png')
    plt.show()

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
    # Check
    print('Value should be: [0.65625, 0.0, 0.65625, 0.65625, '
          '0.65625, -1.3125, 0.79891, -1.45516, 0.79891, -1.3125, '
          '0.65625, 0.65625, -0.79891, 0.14266, -0.79891]')
    print('Value is: ', grad_U([1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0]))

    # fixed points in our model
    origin = np.array([0, 0, 0])
    z2 = 1
    v1 = np.array([0, 0, z2])

    ig = np.array([origin, v1, [1, 0, 0], [0, 1, 0], [1, 1, 0]])
    # ig = np.array([origin, v1])

    # need to first convert the multidimensional array into a vector array
    m, n = np.shape(ig)
    x = m * n
    ig = np.reshape(ig, x)

    # Chilton advised we should use scipy.optimize over his code
    x_min = opt.minimize(U, ig, method='CG', jac=grad_U)

    print(x_min)
    return


def sa4():
    """
    Use one of the functions in SciPy Optimization to find the global minimum
    of your potential function when n=5 using only the potential function itself(
    not the gradient). You cannot use Nelder-Mead for this.
    """
    # fixed points in our model
    origin = np.array([0, 0, 0])
    z2 = 2
    v1 = np.array([0, 0, z2])

    ig = np.array([origin, v1, [1, 0, 0], [0, 1, 0], [1, 1, 0]])
    # ig = np.array([origin, v1])

    # need to first convert the multidimensional array into a vector array
    m, n = np.shape(ig)
    x = m * n
    ig = np.reshape(ig, x)

    # Chilton advised we should use scipy.optimize over his code
    x_min = opt.minimize(U, ig, method='Powell')

    print(x_min)
    return


def sa5():
    """
    Apply the methods used in (1), (3), and (4) when n=6. Rank the methods
    according to reliability and efficiency.
    """
    # QUESTION
    # What do we start each variable in our nodes at?
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
    at the link provided, so your answers can be readily checked. You should do at
    least n=8.
    http://doye.chem.ox.ac.uk/jon/structures/LJ/tables.150.html
    """
    return
