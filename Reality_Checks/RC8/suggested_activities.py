from poisson import *
import numpy as np
from mymesh import mymesh
import matplotlib.pyplot as plt


def sa1():
    """
    Begin with a fin of dimensions 2×2cm, with 1mm thickness (δ). Assume that 5W
    of power is input along the entire left edge, as if the fin were attached to
    dissipate power from a CPU chip with L=2cm side length. Solve the PDE (8.44)
    with M=N=10 subintervals in the x and y directions (11 points). Use the mesh
    command to plot the resulting heat distribution over the xy- plane. What is
    the maximum temperature of the fin, in ◦C?
    """
    xl = 0
    xr = 2
    yb = 0
    yt = 2

    power = 5  # watts
    length = 2  # cm
    width = .1  # cm
    K = 1.68
    H = 0.005
    m = 10
    n = 10

    # np.set_print_options(precision=8,line_width=140)
    w = poisson(xl, xr, yb, yt, m, n, power, length, width, K, H)
    print('w = ', w)
    xp = np.linspace(xl, xr, m + 1)
    yp = np.linspace(yb, yt, n + 1)
    mymesh(xp, yp, w.T, 'x', 'y', 'w')
    return


def sa2():
    """
    Increase the size of the fin to 4×4cm. Input 5W of power along the interval
    [0, 2] on the left side of the fin, as in the previous step. Plot the resulting
    distribution, and find the maximum temperature. Experiment with increased
    values of M and N (use 10, 40, 100). How much does the solution change?
    """
    xl = 0
    xr = 2
    yb = 0
    yt = 2

    power = 5  # watts
    length = 4  # cm
    width = .1  # cm
    K = 1.68
    H = 0.005
    m = 10
    n = 10

    w = poisson(xl, xr, yb, yt, m, n, power, length, width, K, H)
    xp = np.linspace(xl, xr, m + 1)
    yp = np.linspace(yb, yt, n + 1)
    mymesh(xp, yp, w.T, 'x', 'y', 'w')

    lists = np.array([10, 40, 100])
    max_temp = 0
    for i in lists:
        for j in lists:
            m = i
            n = j

            # np.set_print_options(precision=8,line_width=140)
            w = poisson(xl, xr, yb, yt, m, n, power, length, width, K, H)
            print("max temp: ", max(w[0]))
            if max(w[0]) > max_temp:
                max_temp = max(w[0])
            # xp = np.linspace(xl, xr, m + 1)
            # yp = np.linspace(yb, yt, n + 1)
            # mymesh(xp, yp, w.T, 'x', 'y', 'w')

    print("max overall temp: ", max_temp)
    return


def sa3():
    """
    Find the maximum power that can be dissipated by a 4×4cm fin while keep-ing the
    maximum temperature less than 80◦C. Assume the power input is the same as in
    steps 1 and 2.
    """
    xl = 0
    xr = 2
    yb = 0
    yt = 2

    power = 5.4638  # watts
    length = 4  # cm
    width = .1  # cm
    K = 1.68
    H = 0.005
    m = 100
    n = 100
    w = poisson(xl, xr, yb, yt, m, n, power, length, width, K, H)

    while max(w[0]) <= 80.0:
        power += .0001
        w = poisson(xl, xr, yb, yt, m, n, power, length, width, K, H)
        # if (max(w[0]) > 80):
        #     power -= 0.05
        #     w = poisson(xl, xr, yb, yt, m, n, power, length, width, K, H)
        print('power: ', power)
        print('temp: ', max(w[0]))
    # np.set_print_options(precision=8,line_width=140)
    xp = np.linspace(xl, xr, m + 1)
    yp = np.linspace(yb, yt, n + 1)
    mymesh(xp, yp, w.T, 'x', 'y', 'w')
    return


def sa4():
    """
    Replace the aluminum fin with a copper fin, with thermal conductivity
    K=3.85W/cm◦C. Find the maximum power that can be dissipated by a 4×4 cm fin with
    the 2cm power input as in steps 1 and 2, while keeping the maximum temperature
    below 80◦C.
    """
    xl = 0
    xr = 2
    yb = 0
    yt = 2

    power = 5.67  # watts
    length = 4  # cm
    width = .1  # cm
    K = 3.85
    H = 0.005
    m = 100
    n = 100
    w = poisson(xl, xr, yb, yt, m, n, power, length, width, K, H)

    while max(w[0]) <= 80.0:
        power += .001
        w = poisson(xl, xr, yb, yt, m, n, power, length, width, K, H)
        # if (max(w[0]) > 80):
        #     power -= 0.05
        #     w = poisson(xl, xr, yb, yt, m, n, power, length, width, K, H)
        print('power: ', power)
        print('temp: ', max(w[0]))
    # np.set_print_options(precision=8,line_width=140)
    xp = np.linspace(xl, xr, m + 1)
    yp = np.linspace(yb, yt, n + 1)
    mymesh(xp, yp, w.T, 'x', 'y', 'w')
    return


def sa5():
    """
    Plot the maximum power that can be dissipated in step 4 (keeping maximum
    temperature below 80 degrees) as a function of thermal conductivity,
    for 3.67 ≤ K ≤ 5 W/cm◦C
    """
    xl = 0
    xr = 2
    yb = 0
    yt = 2

    power = 5.67  # watts
    length = 4  # cm
    width = .1  # cm
    K = 3.65
    x = np.array([K])
    H = 0.005
    m = 100
    n = 100
    # xx = [1, 2, 3, 4, 5, 6]
    # yy = [3.4, 4.8, 5.7, 9.7, 3.1, 1.4]
    # plt.plot(xx, yy, 'r-', linewidth=1.0, label='test')
    # plt.show()
    w = poisson(xl, xr, yb, yt, m, n, power, length, width, K, H)
    y = np.array([max(w[0])])
    while K < 5.0:
        K += .05
        x = np.append(x, K)
        w = poisson(xl, xr, yb, yt, m, n, power, length, width, K, H)
        y = np.append(y, max(w[0]))
        print('K: ', K)

    print('len x: ', len(x))
    print('len y: ', len(y))
    plt.plot(x, y, 'r-', linewidth=1.0, label='test')
    plt.show()
    plt.savefig("test_fig.png")
    # add axis labels
    return


def sa6():
    """
    Redo step 4 for a water-cooled fin. Assume that water has a convective heat
    transfer coefficient of H=0.1Wcm2◦C, and that the ambient water temperature
    is maintained at 20◦C
    """
    xl = 0
    xr = 2
    yb = 0
    yt = 2

    power = 74.22  # watts
    length = 4  # cm
    width = .1  # cm
    K = 3.85
    H = 0.1
    m = 100
    n = 100
    w = poisson(xl, xr, yb, yt, m, n, power, length, width, K, H)

    while max(w[0]) <= 80.0:
        power += 0.001
        w = poisson(xl, xr, yb, yt, m, n, power, length, width, K, H)
        print('power: ', power)
        print('temp: ', max(w[0]))
    # np.set_print_options(precision=8,line_width=140)
    xp = np.linspace(xl, xr, m + 1)
    yp = np.linspace(yb, yt, n + 1)
    mymesh(xp, yp, w.T, 'x', 'y', 'w')
    """
    Temp found to be 79.9995 at power 74.223
    """
    return


def sa7():
    """
    Using the 4×4cm fin from step 2, try various positions for the 2cm power
    source along the left edge. What is the best location for the power source?
    How much is the maximum temperature reduced compared to the configuration
    in step 2?
    """
    xl = 0
    xr = 2
    yb = 0
    yt = 2

    power = 5  # watts
    length = 4  # cm
    width = .1  # cm
    K = 1.68
    H = 0.005
    m = 10
    n = 10
    x = 1  # shift increment

    w = poisson(xl, xr, yb, yt, m, n, power, length, width, K, H)
    max_w = max(w.max(axis=1))
    print('max w: ', max_w)

    # finds the max of a 2 dimensional array
    shifted_w = shifting_poisson(xl, xr, yb, yt, m, n, power, length, width, K, H, x)
    max_shifted_w = max(shifted_w.max(axis=1))
    print('max shifted w: ', max_shifted_w)

    return

