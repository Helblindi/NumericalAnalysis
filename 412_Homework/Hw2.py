# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 13:27:17 2018

@author: Madison
"""

import numpy as np
import matplotlib.pyplot as plt

"""
Trapezoid Method for Solving Initial Value Problems
where y(t) is a scalar or a (row) vector
f(t,y) evaluates rhs of differential equation y' = f(t, y)
Input: interval [a,b], initial value y0 (np.array, 1 row), number of steps n,
f (a row vector)
Output: time values t, solution y
 """


def trapV(interval, y0, n, f):
    h = float(interval[1] - interval[0]) / n
    t = np.array([i * h for i in np.arange(n + 1)])
    y = np.array([y0])
    d = y0.size
    for i in np.arange(n):
        y = np.append(y, np.reshape(trapstepV(t[i], y[i, :], h, f), (1, d)), axis=0)
    return t, y


def trapstepV(t, w, h, f):
    return w + (h / 2) * (f(t, w) + f(t + h, w + h * f(t, w)))


###############################################################################
print("Computer Problem 6.3.2bc")
"""
Apply the Trapezoid Method with step sizes h = 0.1 and h = 0.01 to
the initial value problems in Exercise 1.  Plot the approximate 
solutions and the correcct solution on [0,1], and find the global
truncation error at t = 1.  Is the reduction in error for h = 0.01
consistent with the order of the Trapezoid Method?

(b) 
y1' = y1 + y2
y2' = -y1 + y2
y1(0) = 1
y2(0) = 0

(c)
y1' = -y1 - y2
y2' = y1 - y2
y1(0) = 1
y2(0) = 0
"""

print("Initial Value Problem b")


# differential equation
def fb(t, y):
    return np.array([-y[0] - y[1], y[0] - y[1]])


# actual solution
def yb(t):
    return np.array([np.exp(-t) * np.cos(t), np.exp(-t) * np.sin(t)])


a = 0.0
b = 1.0
y0 = np.array([1, 0])
h = 0.1
steps = int((b - a) / h)
et1, ey1 = trapV([a, b], y0, steps, fb)

# try step size h = 0.01
h = 0.01
steps = int((b - a) / h)
et2, ey2 = trapV([a, b], y0, steps, fb)

# print global truncation errors
print("Global Truncation Error for (h = 0.1) = [", np.abs(yb(1)[0] - ey1[-1, 0]), np.abs(yb(1)[1] - ey1[-1, 1]), "]")
print("Global Truncation Error for (h = 0.01) = [", np.abs(yb(1)[0] - ey2[-1, 0]), np.abs(yb(1)[1] - ey2[-1, 1]), "]")

fig1 = plt.figure(figsize=(8, 6))
ax1 = fig1.add_subplot(111)
pt = np.linspace(a, b, 50, endpoint=True)
ax1.plot(pt, yb(pt)[0, :], 'b-', label='Exact$y_1$')
ax1.plot(et1, ey1[:, 0], 'ro', label='Approx$y_1$,$h=0.1$')
ax1.plot(et2, ey2[:, 0], 'g--', label='Approx$y_1$,$h=0.01$', linewidth=2.2)
ax1.plot(pt, yb(pt)[1, :], 'r-', label='Exact$y_2$', linewidth=2.2)
ax1.plot(et1, ey1[:, 1], 'go', label='Approx$y_2$,$h=0.1$')
ax1.plot(et2, ey2[:, 1], 'b--', label='Approx$y_2$,$h=0.01$', linewidth=2.2)
plt.xlabel('$t$', fontsize=16)
plt.title('$y_1(t) ,y_2(t)$', fontsize=16)
ax1.legend(fontsize='large', loc='upper right')
plt.show()

print("Initial Value Problem c")


# differential equation
def fb(t, y):
    return np.array([-y[1], y[0]])


# actual solution
def yb(t):
    return np.array([np.cos(t), np.sin(t)])


a = 0.0
b = 1.0
y0 = np.array([1, 0])
h = 0.1
steps = int((b - a) / h)
et1, ey1 = trapV([a, b], y0, steps, fb)

# try step size h = 0.01
h = 0.01
steps = int((b - a) / h)
et2, ey2 = trapV([a, b], y0, steps, fb)

# print global truncation errors
print("Global Truncation Error for (h = 0.1) = [", np.abs(yb(1)[0] - ey1[-1, 0]), np.abs(yb(1)[1] - ey1[-1, 1]), "]")
print("Global Truncation Error for (h = 0.01) = [", np.abs(yb(1)[0] - ey2[-1, 0]), np.abs(yb(1)[1] - ey2[-1, 1]), "]")

fig1 = plt.figure(figsize=(8, 6))
ax1 = fig1.add_subplot(111)
pt = np.linspace(a, b, 50, endpoint=True)
ax1.plot(pt, yb(pt)[0, :], 'b-', label='Exact$y_1$')
ax1.plot(et1, ey1[:, 0], 'ro', label='Approx$y_1$,$h=0.1$')
ax1.plot(et2, ey2[:, 0], 'g--', label='Approx$y_1$,$h=0.01$', linewidth=2.2)
ax1.plot(pt, yb(pt)[1, :], 'r-', label='Exact$y_2$', linewidth=2.2)
ax1.plot(et1, ey1[:, 1], 'go', label='Approx$y_2$,$h=0.1$')
ax1.plot(et2, ey2[:, 1], 'b--', label='Approx$y_2$,$h=0.01$', linewidth=2.2)
plt.xlabel('$t$', fontsize=16)
plt.title('$y_1(t) ,y_2(t)$', fontsize=16)
ax1.legend(fontsize='large', loc='lower right')
plt.show()

###############################################################################
print("Computer Problem 6.3.9")
"""
Adapt orbit.m to solve the two-body problem.  Set the masses to m1 = 0.2, 
m2 = 0.03, and plto the trajectories with initial conditions (x1, y1) = (2,2), 
(x1', y1') = (0.2, -0.2) and (x2, y2) = (0,0), (x2', y2') = (-0.01, 0.01).
"""


# trap step with params
def stepV(t, w, h, params, f):
    return w + (h / 2) * (f(t, w, params) + f(t + h, w + h * f(t, w, params), params))


# equations of motion for the two-body problem
def f(t, x, params):
    g = params[0]
    m1 = params[1]
    m2 = params[2]
    mg1 = m1 * g
    mg2 = m2 * g
    dist = np.sqrt((x[4] - x[0]) ** 2 + (x[6] - x[2]) ** 2)
    z = np.array([x[1], mg2 * (x[4] - x[0]) / dist ** 3, x[3], mg2 * (x[6] - x[2]) / dist ** 3,
                  x[5], mg1 * (x[0] - x[4]) / dist ** 3, x[7], mg1 * (x[2] - x[6]) / dist ** 3])
    return z


def orbit(params, inter, ic, n, p, f):
    # use n points
    h = (inter[1] - inter[0]) / float(n)
    # build t vector
    t = np.zeros((p + 1))
    t[0] = inter[0]
    rt = np.zeros(int(n / p) + 1)
    rt[0] = inter[0]
    # build y vector
    y = np.zeros((p + 1, len(ic)))
    y[0, :] = ic
    ry = np.zeros((int(n / p) + 1, len(ic)))
    ry[0] = ic
    for k in range(1, int(n / p) + 1):
        for i in range(p):
            t[i + 1] = t[i] + h
            y[i + 1, :] = stepV(t[i], y[i, :], h, params, f)
            rt[k] = t[i + 1]
            ry[k, :] = y[i + 1, :]
            t[0] = rt[k]
            y[0:k] = ry[k, :]
    # rt, ry are every pth point in the solution (for plotting)
    return rt, ry


t0 = 0;
tf = 85
x10 = 2;
x1p0 = 0.2;
y10 = 2;
y1p0 = -0.2
x20 = 0.0;
x2p0 = -0.01;
y20 = 0.0;
y2p0 = 0.01

steps = 10000;
pn = 50
g = 1;
m1 = 0.3;
m2 = 0.03

rt, ry = orbit([g, m1, m2], [t0, tf], [x10, x1p0, y10, y1p0, x20, x2p0, y20,
                                       y2p0], steps, pn, f)

fig1 = plt.figure(figsize=(8, 6))
ax1 = fig1.add_subplot(111)
ax1.plot(x10, y10, 'bo')
ax1.plot(x20, y20, 'ro')
ax1.plot(ry[:, 0], ry[:, 2], 'b-', label='Mass 1')
ax1.plot(ry[:, 4], ry[:, 6], 'r-', label='Mass 2')
plt.axes().set_aspect('equal')
ax1.axis([-1.0, 18.0, -15.0, 2.5])
plt.legend(fontsize='large', loc='upper right')
plt.show()

# find ics that produce a circular orbit
t0 = 0;
tf = 85
x10 = 0.0;
x1p0 = 0.0;
y10 = 0.0;
y1p0 = -0.012
x20 = 3.0;
x2p0 = 0.0;
y20 = 0.0;
y2p0 = 0.816
steps = 10000;
pn = 50
g = 1;
m1 = 2.0;
m2 = 0.03

rt, ry = orbit([g, m1, m2], [t0, tf], [x10, x1p0, y10, y1p0, x20, x2p0, y20,
                                       y2p0], steps, pn, f)

fig2 = plt.figure(figsize=(8, 6))
ax2 = fig2.add_subplot(111)
ax2.plot(x10, y10, 'bo')
ax2.plot(x20, y20, 'ro')
ax2.plot(ry[:, 0], ry[:, 2], 'b-', label='Mass 1')
ax2.plot(ry[:, 4], ry[:, 6], 'r-', label='Mass 2')
plt.axes().set_aspect('equal')
ax1.axis([-4.0, 4.0, -4.0, 4.0])
plt.legend(fontsize='large', loc='upper right')
plt.show()


