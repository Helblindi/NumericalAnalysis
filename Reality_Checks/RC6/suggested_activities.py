from tacoma import *
from bisect import *
import matplotlib.pylab as plt
from decimal import Decimal


def sa1():
    print("Suggested Activity 1")
    """
    Run tacoma.py using the trapezoid method with wind speed W=80 km/hr
    and initial conditions y=y′=θ′=0,  θ=0.001. The bridge is stable in
    the torsional dimension if small disturbances in θ die out; unstable
    if they grow far beyond original size. Whicho ccurs for this value of W?
    """
    a = 0
    b = 1000
    n = 50000
    p = 4
    theta0 = 1.0 ** (-3)
    d = 0.01
    W = 80.0
    ot, oy = tacoma_trap([a, b], [0, 0, theta0, 0], n, p, d, W)

    # oy[0] = y, 0y[2] = theta
    # plot for y
    # fig = plt.figure(1, figsize=(8, 6))
    # plt.plot(ot, oy[:, 0], 'b-')
    # plt.xlabel('$t$', fontsize=18, color='Blue')
    # plt.ylabel('$y(t)$', fontsize=18, color='Blue')
    # plt.title('Vertical Displacement')
    # plt.savefig("SA_1A")
    # plt.show()

    # Plot the Torsional Displacement
    fig = plt.figure(2, figsize=(8, 6))
    plt.plot(ot, oy[:, 2], 'b-')
    plt.xlabel('$t$', fontsize=18, color='Blue')
    plt.ylabel('$\\theta(t)$', fontsize=18, color='Blue')
    plt.title('Torsional Displacement')
    plt.savefig("SA_1B")
    plt.show()

    """
    In response to the questions we see that the small disturbances
    in theta die out and thus the bridge is stable in the torsional
    dimension.
    """
    return


def sa2():
    print("Suggested Activity 2")
    """
    Replace the trapezoid method by fourth-order Runge–Kutta to
    improve accuracy. Also, plot y(t) and θ(t), for 0 ≤ t ≤ 1000.
    """
    # Inputs: inter = time interval, ic = [y, y', theta, theta'],
    # n = number of steps, p = steps per point plotted
    a = 0
    b = 1000
    n = 50000
    p = 4
    theta0 = 1.0 ** (-3)
    d = 0.01
    W = 80.0
    ot, oy = tacoma_RK4([a, b], [0, 0, theta0, 0], n, p, d, W)

    # oy[0] = y, 0y[2] = theta
    # plot y
    fig = plt.figure(1, figsize=(8, 6))
    plt.plot(ot, oy[:, 0], 'b-')
    plt.xlabel('$t$', fontsize=18, color='Blue')
    plt.ylabel('$y(t)$', fontsize=18, color='Blue')
    plt.title('Vertical Displacement')
    plt.savefig("SA_2A")
    plt.show()

    # plot theta
    fig = plt.figure(2, figsize=(8, 6))
    plt.plot(ot, oy[:, 2], 'b-')
    plt.xlabel('$t$', fontsize=18, color='Blue')
    plt.ylabel('$\\theta(t)$', fontsize=18, color='Blue')
    plt.title('Torsional Displacement')
    plt.savefig("SA_2B")
    plt.show()

    """
    We see that merely changing the step method did not change the fact
    that the bridge is stable in the torsional dimension.
    """

    return


def sa3():
    print("Suggested Activity 3")
    """
    The system is torsionally stable for W = 50 km/hr. Find the
    magnification factor for a small initial angle. That is, set
    θ(0) = 10-3 and find the ratio of the maximum angle θ(t),
    0 ≤ t < ∞, to θ(0). Is the magnification factor approximately
    consistent for initial angles θ(0)=10-3, 10-4, 10-5, 10-6
    """
    a = 0
    b = 1000
    n = 50000
    p = 4
    theta0 = 10.0 ** (-3)
    d = 0.01
    W = 50.0
    ot, oy = tacoma_RK4([a, b], [0, 0, theta0, 0], n, p, d, W)

    # calculate and print magnification factors for 10**-3, 10**-4, 10**-5, 10**-6
    print('Magnification Factor for theta for theta(0)= %.2E: %2f' % (Decimal(theta0), max(oy[:, 2]) / theta0))

    theta0 = 10.0 ** (-4)
    ot, oy = tacoma_RK4([a, b], [0, 0, theta0, 0], n, p, d, W)
    print('Magnification Factor for theta for theta(0)= %.2E: %2f' % (Decimal(theta0), max(oy[:, 2]) / theta0))

    theta0 = 10.0 ** (-5)
    ot, oy = tacoma_RK4([a, b], [0, 0, theta0, 0], n, p, d, W)
    print('Magnification Factor for theta for theta(0)= %.2E: %2f' % (Decimal(theta0), max(oy[:, 2]) / theta0))

    theta0 = 10.0 ** (-6)
    ot, oy = tacoma_RK4([a, b], [0, 0, theta0, 0], n, p, d, W)
    print('Magnification Factor for theta for theta(0)= %.2E: %2f' % (Decimal(theta0), max(oy[:, 2]) / theta0))

    return


def sa4():
    print("Suggested Activity 4")
    """
    Find the minimum wind speed W for which a small disturbance
    θ(0) = 10 - 3 has a magnification factor of 100 or more. Can
    a consistent magnification factor be defined for this W?
    """
    a = 0
    b = 1000
    n = 50000
    p = 4
    theta0 = 10.0 ** (-3)
    d = 0.01
    W = 59.0
    ot, oy = tacoma_RK4([a, b], [0, 0, theta0, 0], n, p, d, W)

    # calculate and print magnification factors for 10**-3, 10**-4, 10**-5, 10**-6
    print('Magnification Factor for theta for theta(0)= %.2E: %2f' % (Decimal(theta0), max(oy[:, 2]) / theta0))

    # check to see if a consistent magnification factor can be found for this W
    theta0 = 10.0 ** (-4)
    ot, oy = tacoma_RK4([a, b], [0, 0, theta0, 0], n, p, d, W)
    print('Magnification Factor for theta for theta(0)= %.2E: %2f' % (Decimal(theta0), max(oy[:, 2]) / theta0))

    theta0 = 10.0 ** (-5)
    ot, oy = tacoma_RK4([a, b], [0, 0, theta0, 0], n, p, d, W)
    print('Magnification Factor for theta for theta(0)= %.2E: %2f' % (Decimal(theta0), max(oy[:, 2]) / theta0))

    theta0 = 10.0 ** (-6)
    ot, oy = tacoma_RK4([a, b], [0, 0, theta0, 0], n, p, d, W)
    print('Magnification Factor for theta for theta(0)= %.2E: %2f' % (Decimal(theta0), max(oy[:, 2]) / theta0))

    """
    We found W = 59.0 to be the minimum wind speed for which a small
    disturbance θ(0) = 10 - 3 has a magnification factor of 100 or 
    more.  A constant magnification error is found to be ~75.227.
    """
    return


def sa5():
    print("Suggested Activity 5")
    """
    Design and implement a method for computing the minimum wind speed
    in Step 4, to within 0.5 × 10 - 3 km/hr. You may want to use an
    equation solver from Chapter 1.
    """

    def F(W):
        # Setup tacoma bridge with correct initial conditions and interval
        a = 0
        b = 1000
        n = 50000
        p = 4
        theta0 = 10.0 ** (-3)
        d = 0.01
        ot, oy = tacoma_RK4([a, b], [0, 0, theta0, 0], n, p, d, W)

        # Calculate magnification factor and subtract 100 to set up root
        return (max(oy[:, 2]) / theta0) - 100.0

    # a = 58.9 and b = 59.0 bracket the root. Start with those values
    tol = 0.0005
    midpoint, n = bisect(F, 58.9, 59.0, tol)

    print("Found W to be %.4f in %d iterations." % (midpoint, n))

    """
    We used the bisection method to calculate W to be 58.9926 in 7 iterations.
    """

    return


def sa6():
    print("Suggested Activity 6")
    """
    What is the effect of increasing the damping coefficient? Double the current
    value of d (to d = 0.02) and change ω to 3 to adjust for the new d. Compute
    the new critical W (when the magnification factor exceeds 100) and compare
    to the critical W associated with d = 0.01. Can you suggest possible changes
    in design that might have made the bridge less susceptible to torsion?
    """

    def F(W):
        # Setup tacoma bridge with correct initial conditions and interval
        a = 0
        b = 1000
        n = 50000
        p = 4
        theta0 = 10.0 ** (-3)
        d = 0.02
        omega = 3
        ot, oy = tacoma_RK4([a, b], [0, 0, theta0, 0], n, p, d, W, omega)

        # Calculate magnification factor and subtract 100 to set up root
        return (max(oy[:, 2]) / theta0) - 100.0

    tol = 0.0005
    # We ran the bisection method with a = 0.0 and b = 100.0 at first to find the actual root,
    # then we changed the a and b so as to reduce iterations and computation speed.
    # Our original values were a = 0.0, b = 100.0, W = 26.5606 in 17 iterations.
    midpoint, n = bisect(F, 25.0, 27.0, tol)

    print("Found W to be %.4f in %d iterations." % (midpoint, n))

    """
    We used the bisection method to calculate W to be 26.5606 in 17 iterations.
    """

    return
