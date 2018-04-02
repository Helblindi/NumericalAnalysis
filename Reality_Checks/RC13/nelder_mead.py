import numpy as np
import copy


# Program 13.3 Nelder-Mead Search (p. 573)
# Input: function f, best guess xbar (column vector),
# initial search radius rad, and number of steps k
# Output: matrix x whose columns are vertices of simplex
def neldermead(f, xbar, rad, k):
    # in case input array is list or has int data type
    xbar = np.array(xbar, dtype=float)
    # print(xbar)
    n = xbar.shape[0]
    x = np.empty((n, n + 1))
    # each column of x is a simplex vertex
    x[:, 0] = xbar
    x[:, 1:n + 1] = xbar * np.ones((1, n)) + rad * np.eye(n, n)
    # print(x)
    y = np.empty(n + 1)
    for j in range(n + 1):
        # evaluate obj function f at each vertex
        y[j] = f(x[:, j])
    # sort the function values in ascending order
    oy = np.argsort(y)
    y = y[oy]
    # and rank the vertices the same way
    x = x[:, oy]
    #    print(y)
    #    print(x)
    for i in range(k):
        # print('i =', i)
        # xbar is the centroid of the face
        # xbar = np.mean( x[:,0:n], axis=1)
        xbar = np.mean(x[:, 0:n], axis=1)
        # print('x = ', x)
        # omitting the worst vertex xh
        # xh = x[:,n].copy()
        xh = x[:, n]
        xr = 2 * xbar - xh
        yr = f(xr)
        # print(x)
        if yr < y[n - 1]:
            # try expansion xe
            if yr < y[0]:
                xe = 3 * xbar - 2 * xh
                ye = f(xe)
                # accept expansion
                if ye < yr:
                    x[:, n] = xe
                    y[n] = f(xe)
                # accept reflection
                else:
                    x[:, n] = xr
                    y[n] = f(xr)
            # xr is middle of pack, accept reflection
            else:
                x[:, n] = xr
                y[n] = f(xr)
        # xr is still the worst vertex, contract
        else:
            # try outside contraction xoc
            if yr < y[n]:
                xoc = 1.5 * xbar - 0.5 * xh
                yoc = f(xoc)
                # accept outside contraction
                if yoc < yr:
                    x[:, n] = xoc;
                    y[n] = f(xoc)
                # shrink simplex toward best point
                else:
                    for j in range(1, n + 1):
                        x[:, j] = 0.5 * x[:, 0] + 0.5 * x[:, j]
                        y[j] = f(x[:, j])
            # xr is even worse than the previous worst
            else:
                xic = 0.5 * xbar + 0.5 * xh
                yic = f(xic)
                # accept inside contraction
                if yic < y[n]:
                    x[:, n] = xic
                    y[n] = yic
                # shrink simplex toward best point
                else:
                    for j in range(1, n + 1):
                        x[:, j] = 0.5 * x[:, 0] + 0.5 * x[:, j]
                        y[j] = f(x[:, j])

        # re-sort the function values in ascending order
        # and rank the vertices the same way
        oy = np.argsort(y)
        y = y[oy]
        x = x[:, oy]
    # print(i,y)
    return x
