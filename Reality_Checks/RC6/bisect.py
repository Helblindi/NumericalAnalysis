import numpy as np


# bisect(f, a, b, tol) computes approximate solution of f(x) = 0
# Input: function f, a,b such that f(a)*f(b)<0, tolerance tol
# Output: Approximate solution xc and number of iterations
def bisect(f, a, b, tol):
    # evaluate f(x) at a and b
    fa = f(a)
    fb = f(b)
    n = 0

    # check to ensure root is bracketed by [a, b]
    if np.sign(fa * fb) >= 0:
        print('f(a)f(b) < 0 not satisfied!')
        quit()

    # repeat loop until the interval [a, b] is small enough
    while (b - a) / 2. > tol:

        # count number if iterations
        n = n + 1

        # compute centerpoint of interval [a, b]
        c = (a + b) / 2.

        # evaluate f at centerpoint
        fc = f(c)

        # if c is a soluiton, done
        if fc == 0:
            return c

        # check if a and c bracket root, if yes, then b = c, otherwise a = c
        # the new interval is [a, c], or b = c
        if np.sign(fc * fa) < 0:
            b = c
            fb = fc

        # the new interval is [c, b], or a = c
        else:
            a = c
            fa = fc

    # return approximate root and number of iterations
    # midpoint is best approximation
    return [(a + b) / 2., n]


# run bisect() on an example
# define a test function
def f(x):
    return np.cos(x) - x


def test():
    # define error tolerance
    tol = 1.0e-8

    # run bisect()
    res = bisect(f, 0, 1, tol)

    # print results, the first thing returned is the approximate root,
    # the second thing is the number of iterations
    print('approximate root = %14.12f' % res[0])
    print('number of iterations = %d' % res[1])