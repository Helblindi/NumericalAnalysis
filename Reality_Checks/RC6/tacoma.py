# tacoma.py
# implements the code for RC6
# Inputs: inter = time interval, ic = [y, y', theta, theta'],
# n = number of steps, p = steps per point plotted
# Calls a one-step method such as trapstep or RK4step

import numpy as np


def tacoma_RK4(inter, ic, n, p, d, W, omega=2 * np.pi * 38 / 60):
    # use n points
    h = (inter[1] - inter[0]) / float(n)

    # build t vectors
    t = np.zeros((p + 1))
    t[0] = inter[0]
    rt = np.zeros(int(n / p) + 1)
    rt[0] = inter[0]

    # build y vector
    # y[0] is vertical displacement
    # y[1] is the derivative of y[0]
    # y[2] is the rotation of the roadbed
    # y[3] is the derivative of y[2]
    y = np.zeros((p + 1, 4))
    y[0, :] = ic
    ry = np.zeros((int(n / p) + 1, 4))
    ry[0] = ic

    # compute solution
    for k in np.arange(1, int(n / p) + 1):
        for i in np.arange(p):
            t[i + 1] = t[i] + h
            y[i + 1, :] = RK4step(t[i], y[i, :], h, d, W, omega)
        rt[k] = t[i + 1]
        ry[k, :] = y[i + 1, :]
        t[0] = rt[k]
        y[0:k] = ry[k, :]

    # return time and solution
    return rt, ry


def tacoma_trap(inter, ic, n, p, d, W, omega=2 * np.pi * 38 / 60):
    # use n points
    h = (inter[1] - inter[0]) / float(n)

    # build t vectors
    t = np.zeros((p + 1))
    t[0] = inter[0]
    rt = np.zeros(int(n / p) + 1)
    rt[0] = inter[0]

    # build y vector
    # y[0] is vertical displacement
    # y[1] is the derivative of y[0]
    # y[2] is the rotation of the roadbed
    # y[3] is the derivative of y[2]
    y = np.zeros((p + 1, 4))
    y[0, :] = ic
    ry = np.zeros((int(n / p) + 1, 4))
    ry[0] = ic

    # compute solution
    for k in np.arange(1, int(n / p) + 1):
        for i in np.arange(p):
            t[i + 1] = t[i] + h
            y[i + 1, :] = RK4step(t[i], y[i, :], h, d, W, omega)
        rt[k] = t[i + 1]
        ry[k, :] = y[i + 1, :]
        t[0] = rt[k]
        y[0:k] = ry[k, :]

    # return time and solution
    return rt, ry


# trapezoid method step
def trapstep(t, w, h, d, W, omega):
    return w + (h / 2) * (ydot(t, w, d, W, omega) + ydot(t + h, w + h * ydot(t, w, d, W, omega), d, W, omega))


# Runge Kutta step
def RK4step(t, w, h, d, W, omega):
    s1 = ydot(t, w, d, W, omega)
    s2 = ydot(t + h / 2, w + (h / 2) * s1, d, W, omega)
    s3 = ydot(t + h / 2, w + (h / 2) * s2, d, W, omega)
    s4 = ydot(t + h, w + h * s3, d, W, omega)
    return w + (h / 6) * (s1 + 2 * s2 + 2 * s3 + s4)


# this is f(t,y) in y' = f(t, y)
def ydot(t, y, d, W, omega):
    K = 1000
    m = 2500
    L = 6
    a = 0.2
    # omega = 2 * np.pi * 38 / 60
    c1 = K / (m * a)
    a1 = np.exp(a * (y[0] - L * np.sin(y[2])))
    a2 = np.exp(a * (y[0] + L * np.sin(y[2])))
    ydot = np.zeros(4)
    ydot[0] = y[1]
    ydot[1] = -d * y[1] - c1 * (a1 + a2 - 2) + 0.2 * W * np.sin(omega * t)
    ydot[2] = y[3]
    ydot[3] = -d * y[3] + c1 * (3 / L) * np.cos(y[2]) * (a1 - a2)
    return ydot


# a = 0
# b = 1000
# n = 50000
# p = 4
# theta0 = 10.0 ** (-3)
# ot, oy = tacoma([a, b], [0, 0, theta0, 0], n, p)
#
# # oy[0] = y, 0y[2] = theta
# fig = plt.figure(1, figsize=(8, 6))
# plt.plot(ot, oy[:, 0], 'b-')
# plt.xlabel('$t$', fontsize=18, color='Blue')
# plt.ylabel('$y(t)$', fontsize=18, color='Blue')
# plt.title('Vertical Displacement')
#
# fig = plt.figure(2, figsize=(8, 6))
# plt.plot(ot, oy[:, 2], 'b-')
# plt.xlabel('$t$', fontsize=18, color='Blue')
# plt.ylabel('$\\theta(t)$', fontsize=18, color='Blue')
# plt.title('Torsional Displacement')
# plt.show()
