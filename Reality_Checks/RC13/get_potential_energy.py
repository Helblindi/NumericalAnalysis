import numpy as np


def U(points):
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


def test():
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
