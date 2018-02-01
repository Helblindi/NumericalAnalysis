# -*- coding: utf-8 -*-
"""
Created on Mon Jan 29 11:05:33 2018

@author: Madison
"""

import numpy as np
import matplotlib.pyplot as plt

print("Computer Problem 6.5.1f")
"""
Write a python implementation of RK23 (Example 6.19), and apply to
approximating the solutions of the IVPs in Exercise 6.1.3 with a 
relative tolerance of 10^-8 on [0,1].  Ask the program to stop exactly
at the endpoint t=1.  Report the maximum step size used and the number
of steps.
"""
# Solve example 6.1 (p. 283) using RK23 (example 6.19) with a relative
# tolerance of 10^(-8) on [0, 1]. Stop exactly at t = 1.
# Report maximum step size used and the number of steps.

# Example 6.1
def f(t, y):
    return t * y + t ** 3


def y(t):
    return 3 * np.exp((t ** 2) / 2) - t ** 2 - 2


# One step of Scalar RK23 Method
# Input: current time t, current value w, stepsize h, RHS f
# Output: new h, and approximate solution w at time t + h
def RK23stepVf(t, w, h, f):
    tol = 10 ** (-5.0)

    # Get each piece needed for the RK4 method
    s1 = f(t, w)
    s2 = f(t + h, w + h * s1)
    s3 = f(t + h / 2, w + (h / 4) * (s1 + s2))
    error = (h / 3) * np.abs(s1 - 2 * s3 + s2)

    # compare error to relative tolerance
    while (error > tol * np.abs(w)):
        # recalculate our necessary values
        h = h / 2
        s1 = f(t, w)
        s2 = f(t + h, w + h * s1)
        s3 = f(t + h / 2, w + (h / 4) * (s1 + s2))
        error = (h / 3) * np.abs(s1 - 2 * s3 + s2)

    w = w + (h / 6) * (s1 + 4 * s2 + s3)

    h0 = h
    if (error < tol * np.abs(w) / 10):
        h0 = 2 * h;
    return h, w, h0


a = 0
b = 1
y0 = 1

h0 = (b - a) / 2.0

t = np.array([a])
w = np.array([y0])
count = 0

while t[count] < b:
    count = count + 1

    nh, nw, h0 = RK23stepVf(t[count - 1], w[count - 1], h0, f)

    w = np.append(w, nw)
    t = np.append(t, t[count - 1] + nh)

if t[count] > b:
    h = b - t[count - 1]
    nh, nw, n0 = RK23stepVf(t[count - 1], w[count - 1], h0, f);
    w = np.append(w, nw)
    t = np.append(t, t[count - 1] + nh)

# display the meximum step size and the number of steps
print('max step size = %f' % np.max(np.diff(t)))
print('number of steps = %d' % (t.size - 1))

# plot our results
pt = np.linspace(a, b, 200, endpoint=True)
plt.plot(pt, y(pt), 'b-')
plt.plot(t, w, 'go-')
plt.plot(t, 0.95 * np.ones(t.size), 'ro', markersize=3.0)
plt.ylim([0.9, 2.0])

plt.show()

###############################################################################
print("Computer Problem 6.6.1b")
"""
Apply Backward Euler, using Newton's Method as a solver, for the intial value
problems.  Which of the equlibrium solutions are approached by the approximate
solution? Apply Euler's Method.  For what approximate rande of h can Euler be 
used successfully to converge to the equilibrium?  Plot approximate solutions 
given by Backward Euler, and by Euler with an excessive step size.
y' = 6y - 6y^2
y(0) = 1/2
t in [0,20]
"""