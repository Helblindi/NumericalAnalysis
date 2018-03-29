import numpy as np


# def get_distance(r_i, r_j):
#     return np.sqrt((r_i[0] - r_j[0])**2 +
#                    (r_i[1] - r_j[1])**2 +
#                    (r_i[2] - r_j[2])**2)
def get_distance(x, y):
    return np.sqrt((x-y)**2)


def get_potential_energy(points):
    potential_energy = 0
    for i in range(0, len(points)):
        for j in range(i + 1, len(points)):
            # print('i: ', i)
            # print('j: ', j)
            # print('i+j: ', i+j)
            r_ij = get_distance(points[i], points[j])
            potential_energy += 1 / (r_ij**12) - 2 / (r_ij**6)
    return potential_energy


def test():
    two_points = np.array([[0, 0, 0], [0, 0, 1]])
    z2 = 1
    origin = np.array([0, 0, 0])
    v1 = np.array([0, 0, z2])
    points_array = np.array([origin, v1, [3, 3, 3], [4, 4, 4], [5, 5, 5]])
    print('Potential Energy: ', get_potential_energy(two_points))
    return
