import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def U(points):
    """
    Calculates the potential energy of a given orientation of nodes
    :param points: single dimensional array representing a node configuration
    :return: potential energy of the given node configuration
    """
    # Calculates potential energy of given orientation of nodes
    # U(1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0) = -6.0102023319615911
    # calculate how many nodes we have
    n = int(len(points) / 3)

    # convert the points vector back into multidimensional array
    points = np.reshape(points, (n, 3))

    # calculate potential energy
    potential_energy = 0
    for i in range(0, n - 1):
        for j in range(i + 1, n):
            # print('i: ', i)
            # print('j: ', j)
            # print('i+j: ', i+j)
            r_ij = np.linalg.norm(points[i] - points[j])
            if r_ij != 0:
                potential_energy += 1 / (r_ij**12) - 2 / (r_ij**6)

                # is this way to increment potential energy more efficient?
                # potential_energy += (1 / (r_ij**6)) * (1 / (r_ij ** 6) - 2)

    # return potential energy
    return potential_energy


def grad_U(points):
    """
    Gradient functions for U
    :param points: single dimensional array representing a node configuration
    :return: gradient corresponding to the node configuration
    """
    """
    Check:
    ∇U(1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0) = [0.65625, 0.0,
    0.65625, 0.65625,0.65625,-1.3125, 0.79891, -1.45516, 0.79891,-1.3125,
    0.65625, 0.65625,-0.79891, 0.14266,-0.79891]
    """

    # convert the points vector back into multidimensional array
    n = int(len(points) / 3)
    points = np.reshape(points, (n, 3))

    # fill gradient array
    gradient = np.zeros((n, 3))
    for i in range(0, n - 1):
        for j in range(i + 1, n):
            r_ij = np.linalg.norm(points[i] - points[j])

            if r_ij > 0:
                # derivative of the Lennard-Jones potential
                d_ljp = -12 / (r_ij**13) + 12 / (r_ij**7)

                # fill gradient matrix
                for it in range(0, 3):
                    # print('i: ', i, ' j: ', j, ' it: ', it)
                    gradient[i][it] += d_ljp * (points[i][it] - points[j][it]) / r_ij
                    gradient[j][it] += d_ljp * (points[j][it] - points[i][it]) / r_ij

    # need to first convert the multidimensional array into a vector array
    m, n = np.shape(gradient)
    x = m * n
    gradient = np.reshape(gradient, x)

    # return gradient
    return gradient


def translate_and_rotate(points):
    """
    Translates and rotates a node configuration in order to ensure that the
    first point is the origin, and that the second points is fixed on the z-axis.
    :param points: single dimensional array representing a node configuration
    :return: points - a translated and rotated node configuration that is in
    line with Reality Check 13
    """
    # convert the points vector back into multidimensional array
    n = int(len(points) / 3)
    points = np.reshape(points, (n, 3))

    # Transform all points so first point is at the origin
    origin = points[0]
    for i in range(1, n - 1):
        points[i] = points[i] - origin

    # Set the first point to the origin
    points[0] = np.array([0, 0, 0])

    # Now rotate the points in y direction to as to set x value of the second point to 0
    # x = x*cos(t) - z*sin(t)
    # z = z*cos(t) + x*sin(t)
    x = points[1][0]
    z = points[1][2]

    # first calculate angle of rotation in radians
    theta = np.arctan(x / z)

    cos = np.cos(theta)
    sin = np.sin(theta)

    # now create the matrix of rotation
    A = np.array([[cos, 0, -sin], [0, 1, 0], [sin, 0, cos]])

    # apply rotation matrix to each point
    for i in range(0, len(points)):
        points[i] = np.dot(A, points[i])

    # account for machine error
    if points[1][0] < .0005:
        points[1][0] = 0.

    # Now rotate the points in x direction to as to set y value of the second point to 0
    # y = y*cos(t) - z*sin(t)
    # z = z*cos(t) + y*sin(t)
    y = points[1][1]
    z = points[1][2]

    # first calculate angle of rotation in radians
    theta = np.arctan(y / z)

    cos = np.cos(theta)
    sin = np.sin(theta)

    # now create the matrix of rotation
    A = np.array([[1, 0, 0], [0, cos, -sin], [0, sin, cos]])

    # apply rotation matrix to each point
    for i in range(0, len(points)):
        points[i] = np.dot(A, points[i])

    # account for machine error
    if points[1][1] < .0005:
        points[1][1] = 0.

    # lastly, convert the multidimensional array into a vector array
    m, n = np.shape(points)
    x = m * n
    points = np.reshape(points, x)

    return points


def plot_configuration(points, plot_name='plot', plot_title='plot_title'):
    """
    Plots node configuration
    :param points: single dimensional array of given node configuration
    :param plot_name: name plot to be saved as
    :param plot_title: title of plot
    :return: saved plot
    """
    # first to convert points to a multidimensional array
    n = int(len(points) / 3)
    print('n: ', n)
    points = np.reshape(points, (n, 3))

    # now create empty array for each coordinate
    x = np.array([])
    y = np.array([])
    z = np.array([])

    # loop through the points and store each coordinate
    # value into its corresponding coordinate array
    for i in range(0, n):
        x = np.append(x, [points[i][0]])
        y = np.append(y, points[i][1])
        z = np.append(z, points[i][2])

    fig = plt.figure()
    ax = Axes3D(fig)

    # plot each point
    ax.scatter(x, y, z, 'bo')

    # plot line connecting the points
    # https://stackoverflow.com/questions/41563347/connecting-points-to-a-central-point-on-3d-scatter-python?rq=1
    for a, b, c in zip(x, y, z):
        for d, e, f in zip(x, y, z):
            ax.plot3D([a, d], [b, e], [c, f], 'b-')

    plt.title(plot_title, fontsize=24)
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')
    plt.savefig(plot_name)
    plt.show()


def test():
    """
    function to test U
    :return:
    """
    two_points = np.array([[0, 0, 0], [0, 0, 1]])
    z2 = 1
    origin = np.array([0, 0, 0])
    v1 = np.array([0, 0, z2])
    points_array = np.array([origin, v1, [3, 3, 3], [4, 4, 4], [5, 5, 5]])
    print('Potential Energy: ', U(two_points))
    # U(1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0) = -6.0102023319615911
    print(U([1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0]))
    # check out
    return
