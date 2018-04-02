import numpy as np


def get_potential_energy(points):
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


def test():
    two_points = np.array([[0, 0, 0], [0, 0, 1]])
    z2 = 1
    origin = np.array([0, 0, 0])
    v1 = np.array([0, 0, z2])
    points_array = np.array([origin, v1, [3, 3, 3], [4, 4, 4], [5, 5, 5]])
    print('Potential Energy: ', get_potential_energy(two_points))
    # U(1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0) = -6.0102023319615911
    print(get_potential_energy([1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0]))
    # check out
    return
