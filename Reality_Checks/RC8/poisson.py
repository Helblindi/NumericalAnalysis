# Program 8.5 Finite difference solver for 2D Poisson equation
# with Dirichlet boundary conditions on a rectangle

# Example 8.8
# Laplacian u = f(x,y) on [0,1] x [1,2] (X x Y)
# u(x,1) = ln(x^2 + 1)
# u(x,2) = ln(x^2 + 4)
# u(0,y) = 2 ln(y)
# u(1,y) = ln(y^2 + 1)
# exact solution is u(x,y) = ln(x^2 + y^2)

# Input: rectangle domain [xl,xr]x[yb,yt] with MxN space steps and
#        data functions f, g1, g2, g3, g4
# Output: matrix w holding solution values
# Example usage: w = poisson(xl, xr, yb, yt, f, u(x,1), u(x,2),
#                    u(0,y), u(1,y), M, N)

import numpy as np
from mymesh import mymesh


def poisson(xl, xr, yb, yt, M, N, Power, length, width, K, H):
    P = Power  # watts
    L = length  # fin length in cm
    D = width  # find width in mm
    K = K  # w/cm Celsius thermal conductivity
    H = H  # convective heat transfer coefficient w/cm^2 Celsius
    F = 2 * H / float(K * D)
    m = M + 1
    n = N + 1
    mn = m * n
    h = float((xr - xl) / M)
    h2 = h ** 2
    k = float((yt - yb) / N)
    k2 = k ** 2
    # set up mesh points
    x = np.linspace(xl, xr, m)
    y = np.linspace(yb, yt, n)
    # build A and b
    A = np.zeros((mn, mn))
    b = np.zeros(mn)
    # interior points
    for i in range(1, m - 1):
        for j in range(1, n - 1):
            A[i + j * m, i - 1 + j * m] = 1. / h2
            A[i + j * m, i + 1 + j * m] = 1. / h2
            A[i + j * m, i + j * m] = -2. / h2 - 2. / k2 - F
            A[i + j * m, i + (j - 1) * m] = 1. / k2
            A[i + j * m, i + (j + 1) * m] = 1. / k2
            b[i + j * m] = 0

    for i in range(m):
        # bottom
        j = 0
        A[i + j * m, i + j * m] = -3 / (2 * k) + (H / K)
        A[i + j * m, i + (j + 1) * m] = 4 / (2 * k)
        A[i + j * m, i + (j + 2) * m] = -1 / (2 * k)
        b[i + j * m] = 0
        # top
        j = n - 1
        A[i + j * m, i + j * m] = -3 / (2 * k) + (H / K)
        A[i + j * m, i + (j - 1) * m] = 4 / (2 * k)
        A[i + j * m, i + (j - 2) * m] = -1 / (2 * k)
        b[i + j * m] = 0
    # left and right boundaries
    for j in range(1, n - 1):
        # left
        i = 0
        A[i + j * m, i + 1 + j * m] = -3 / (2 * h)
        A[i + j * m, i + 2 + j * m] = 4 / (2 * h)
        A[i + j * m, i + 3 + j * m] = -1 / (2 * h)
        b[i + j * m] = -P/(L*D*K)
        # right
        i = m - 1
        A[i + j * m, i + j * m] = -3 / (2 * h) + (H / K)
        A[i + j * m, i - 1 + j * m] = 4 / (2 * h)
        A[i + j * m, i - 2 + j * m] = -1 / (2 * h)
        b[i + j * m] = 0

    # matprint(A)
    # solve for v
    v = np.linalg.solve(A, b)
    # translate form v to w
    w = np.reshape(v, (m, n), order='F')
    return w


def shifting_poisson(xl, xr, yb, yt, M, N, Power, length, width, K, H, shift):
    # x is shift of power source, where x = 0 represents standard used in other problems
    P = Power  # watts
    L = length  # fin length in cm
    D = width  # find width in mm
    K = K  # w/cm Celsius thermal conductivity
    H = H  # convective heat transfer coefficient w/cm^2 Celsius
    F = 2 * H / float(K * D)
    m = M + 1
    n = N + 1
    mn = m * n
    h = float((xr - xl) / M)
    h2 = h ** 2
    k = float((yt - yb) / N)
    k2 = k ** 2
    # set up mesh points
    x = np.linspace(xl, xr, m)
    y = np.linspace(yb, yt, n)
    # build A and b
    A = np.zeros((mn, mn))
    b = np.zeros(mn)
    # interior points
    for i in range(1, m - 1):
        for j in range(1, n - 1):
            A[i + j * m, i - 1 + j * m] = 1. / h2
            A[i + j * m, i + 1 + j * m] = 1. / h2
            A[i + j * m, i + j * m] = -2. / h2 - 2. / k2 - F
            A[i + j * m, i + (j - 1) * m] = 1. / k2
            A[i + j * m, i + (j + 1) * m] = 1. / k2
            b[i + j * m] = 0

    for i in range(m):
        # bottom
        j = 0
        A[i + j * m, i + j * m] = -3 / (2 * k) + (H / K)
        A[i + j * m, i + (j + 1) * m] = 4 / (2 * k)
        A[i + j * m, i + (j + 2) * m] = -1 / (2 * k)
        b[i + j * m] = 0
        # top
        j = n - 1
        A[i + j * m, i + j * m] = -3 / (2 * k) + (H / K)
        A[i + j * m, i + (j - 1) * m] = 4 / (2 * k)
        A[i + j * m, i + (j - 2) * m] = -1 / (2 * k)
        b[i + j * m] = 0
    # left and right boundaries
    for j in range(1, n - 1):
        # left
        i = 0
        A[i + j * m, i + 1 + shift + j * m] = -3 / (2 * h)
        A[i + j * m, i + 2 + shift + j * m] = 4 / (2 * h)
        A[i + j * m, i + 3 + shift + j * m] = -1 / (2 * h)
        b[i + j * m] = -P/(L*D*K)
        # right
        i = m - 1
        A[i + j * m, i + j * m] = -3 / (2 * h) + (H / K)
        A[i + j * m, i - 1 + j * m] = 4 / (2 * h)
        A[i + j * m, i - 2 + j * m] = -1 / (2 * h)
        b[i + j * m] = 0

    # matprint(A)
    # solve for v
    v = np.linalg.solve(A, b)
    # translate form v to w
    w = np.reshape(v, (m, n), order='F')
    return w


def matprint(mat, fmt="g"):
    col_maxes = [max([len(("{:"+fmt+"}").format(x)) for x in col]) for col in mat.T]
    for x in mat:
        for i, y in enumerate(x):
            print(("{:"+str(col_maxes[i])+fmt+"}").format(y), end="  ")
        print("")
    return


# # define data
# # this is Example 8.8
# # f is the external forcing function
# def f(x, y):
#     return 0
#
#
# # g1 is the boundary condition on the bottom, y = 1
# def g1(x):
#     return np.log(x ** 2 + 1)
#
#
# # g2 is the boundary condition on the top, y = 2
# def g2(x):
#     return np.log(x ** 2 + 4)
#
#
# # g3 is the boundary condition on the left, x = 0
# def g3(y):
#     return 2 * np.log(y)
#
#
# # g4 is the boundary condition on the right, x = 1
# def g4(y):
#     return np.log(y ** 2 + 1)
#
#
# def u_exact(x, y):
#     return np.log(x ** 2 + y ** 2)


# Everything should be in a function
# Sample code that Dr. Chilton gave us
def sample():
    xl = 0
    xr = 1
    yb = 1
    yt = 2

    power = 5  # watts
    m = 20
    n = 20

    # np.set_print_options(precision=8,line_width=140)
    w = poisson(xl, xr, yb, yt, m, n, power)
    # print('w = ', w)
    # xp = np.linspace(xl, xr, m + 1)
    # yp = np.linspace(yb, yt, n + 1)
    # mymesh(xp, yp, w.T, 'x', 'y', 'w')
    #
    # mx, my = np.meshgrid(xp, yp)
    # ex = u_exact(mx, my)
    # error = np.abs(ex - w.T)
    # mymesh(xp, yp, error, 'x', 'y', 'w')
    return 0
