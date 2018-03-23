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


def poisson(xl, xr, yb, yt, M, N, Power, length, width, K, H, p_int=[0,2]):
    P = Power  # watts
    L = length  # fin length in cm
    D = width  # find width in mm
    # K w/cm Celsius thermal conductivity
    # H convective heat transfer coefficient w/cm^2 Celsius
    F = 2 * H / float(K * D)
    m = M + 1
    n = N + 1
    mn = m * n
    h = float((xr - xl) / M)
    h2 = h ** 2
    k = float((yt - yb) / N)
    k2 = k ** 2

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

    for i in range(m):
        # bottom
        j = 0
        A[i + j * m, i + j * m] = -3 / (2 * k) + (H / K)
        A[i + j * m, i + (j + 1) * m] = 4 / (2 * k)
        A[i + j * m, i + (j + 2) * m] = -1 / (2 * k)
        # top
        j = n - 1
        A[i + j * m, i + j * m] = -3 / (2 * k) + (H / K)
        A[i + j * m, i + (j - 1) * m] = 4 / (2 * k)
        A[i + j * m, i + (j - 2) * m] = -1 / (2 * k)

    # left and right boundaries
    for j in range(1, n - 1):
        # left
        i = 0
        A[i + j * m, i + j * m] = -3 / (2 * h) + H / K
        A[i + j * m, i + 1 + j * m] = 4 / (2 * h)
        A[i + j * m, i + 2 + j * m] = -1 / (2 * h)

        if j * h >= p_int[0] and j * h <= p_int[1]:
            A[i + j * m, i + j * m] = -3 / (2 * h)
            b[i + j * m] = -P / ((p_int[1]-p_int[0]) * D * K)
        # right
        i = m - 1
        A[i + j * m, i + j * m] = -3 / (2 * h) + (H / K)
        A[i + j * m, i - 1 + j * m] = 4 / (2 * h)
        A[i + j * m, i - 2 + j * m] = -1 / (2 * h)

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
    return
