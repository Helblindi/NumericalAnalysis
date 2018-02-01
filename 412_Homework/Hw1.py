# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 15:26:53 2018

@author: Madison
"""

import numpy as np
import matplotlib.pyplot as plt

"""
Euler Method for Solving Initial Value Problems
where y(t) is a scalar or a (row) vector
f(t,y) evaluates rhs of differential equation y' = f(t, y)
Input: interval [a,b], initial value y0 (np.array, 1 row), number of steps n,
f (a row vector)
Output: time values t, solution y
 """


def eulerV(interval, y0, n, f):
    h = float(interval[1] - interval[0]) / n
    t = np.array([i * h for i in np.arange(n + 1)])
    y = np.array([y0])
    d = y0.size
    for i in np.arange(n):
        y = np.append(y, np.reshape(eulerstepV(t[i], y[i, :], h, f), (1, d)), axis=0)
    return t, y


def eulerstepV(t, w, h, f):
    return w + h * f(t, w)


###############################################################################
print("6.1 3ab")
"""
Plot the Euler's Method approcimate solutions for the IVPs in Exercise 4 
on [0,1] for step sizes h = 0.1, 0.05, and 0.025, along with the exact solutions
"""
# Variables that will be used for both (a) and (c)
interval = [0, 1]
y0 = np.zeros(1)
n1 = 10
n2 = 20
n3 = 40


# a) y' = t + y; y(0) = 0
def fa(t, y):
    return t + y


def solutionA(t):
    return np.exp(t) - t - 1


ta = np.linspace(0, 1, 100)

# n1
t1, y1 = eulerV(interval, y0, n1, fa)
plt.plot(t1, y1, 'bo')
plt.plot(t1, y1, 'r-')
plt.plot(ta, solutionA(ta))
plt.title('(a) with h = 0.1')
plt.show()

# n2
t1, y1 = eulerV(interval, y0, n2, fa)
plt.plot(t1, y1, 'bo')
plt.plot(t1, y1, 'r-')
plt.plot(ta, solutionA(ta))
plt.title('(a) with h = 0.05')
plt.show()

# n3
t1, y1 = eulerV(interval, y0, n3, fa)
plt.plot(t1, y1, 'bo')
plt.plot(t1, y1, 'r-')
plt.plot(ta, solutionA(ta))
plt.title('(a) with h = 0.025')
plt.show()


# c) y' = 4t - 2y; y(0) = 0
def fc(t, y):
    return 4 * t - 2 * y


def solutionC(t):
    return 2 * t - 1 - np.exp(-2 * t)


tc = np.linspace(0, 1, 100)

# n1
t1, y1 = eulerV(interval, y0, n1, fc)
plt.plot(t1, y1, 'bo')
plt.plot(t1, y1, 'r-')
plt.plot(tc, solutionC(tc))
plt.title('(a) with h = 0.1')
plt.show()

# n2
t1, y1 = eulerV(interval, y0, n2, fc)
plt.plot(t1, y1, 'bo')
plt.plot(t1, y1, 'r-')
plt.plot(tc, solutionC(tc))
plt.title('(a) with h = 0.05')
plt.show()

# n3
t1, y1 = eulerV(interval, y0, n3, fc)
plt.plot(t1, y1, 'bo')
plt.plot(t1, y1, 'r-')
plt.plot(tc, solutionC(tc))
plt.title('(a) with h = 0.025')
plt.show()

###############################################################################
print("6.1 6ac")
"""
For the initial value problems in Exercise 4, make a log-log plot of the error 
of Euler's Method at t = 2 as a function of h = 0.1 x 2^-k for 0<=k<=5 
(also slope)
"""
# a) y' = t + y; y(0) = 0
# Solution: y(t) = e^t - t - 1

k = np.linspace(0, 5, 50)
h = 0.1 * 2 ** (-k)
steps = 2 / h
print(steps)
count = 0
ma = np.zeros(len(steps))
for i in steps:
    t, y = eulerV([0, 2], y0, int(i), fa)
    error = np.abs(solutionA(2) - y[int(i)])
    ma[count] = error
    count += 1

plt.plot(steps, ma, 'bo')
plt.plot(steps, ma, 'r-')
plt.title('log-log plot of error (a)')
plt.show()

# c) y' = 4t - 2y; y(0) = 0
# Solution: y(t) = 2t - 1 - e^-2t

count = 0
mc = np.zeros(len(steps))
for i in steps:
    t, y = eulerV([0, 2], y0, int(i), fc)
    error = np.abs(solutionC(2) - y[int(i)])
    mc[count] = error
    count += 1

plt.plot(steps, mc, 'bo')
plt.plot(steps, mc, 'r-')
plt.title('log-log plot of error (c)')
plt.show()
###############################################################################