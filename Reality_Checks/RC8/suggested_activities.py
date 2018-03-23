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
    u_b = 20.0

    power = 5  # watts
    length = 2  # cm
    width = .1  # cm
    K = 1.68
    H = 0.005
    m = 10
    n = 10

    # np.set_print_options(precision=8,line_width=140)
    w = poisson(xl, xr, yb, yt, m, n, power, length, width, K, H) + u_b
    xp = np.linspace(xl, xr, m + 1)
    yp = np.linspace(yb, yt, n + 1)
    mymesh(xp, yp, w.T, 'x (cm)', 'y (cm)', 'temperature (◦C)', 'SA 1 Heat Distribution')
    max_w = max(w.max(axis=1))
    print('max w: ', max_w)
    """
    We found the maximum temperature across the xy plane given the conditions to be 164.936.
    """
    return


def sa2():
    """
    Increase the size of the fin to 4×4cm. Input 5W of power along the interval
    [0, 2] on the left side of the fin, as in the previous step. Plot the resulting
    distribution, and find the maximum temperature. Experiment with increased
    values of M and N (use 10, 40, 100). How much does the solution change?
    """
    xl = 0
    xr = 4
    yb = 0
    yt = 4
    u_b = 20.0

    power = 5  # watts
    length = 4  # cm
    width = .1  # cm
    K = 1.68
    H = 0.005
    m = 10
    n = 10

    w = poisson(xl, xr, yb, yt, m, n, power, length, width, K, H) + u_b
    xp = np.linspace(xl, xr, m + 1)
    yp = np.linspace(yb, yt, n + 1)
    mymesh(xp, yp, w.T, 'x', 'y', 'w')

    lists = np.array([10, 40, 100])
    max_temp = 100.0
    for i in lists:
        m = i
        n = i

        # np.set_print_options(precision=8,line_width=140)
        w = poisson(xl, xr, yb, yt, m, n, power, length, width, K, H)
        print("m = {}, n = {}, max temp: {}".format(m, n, max(w[0])))
        if max(w[0]) < max_temp:
            max_temp = max(w[0])
        xp = np.linspace(xl, xr, m + 1)
        yp = np.linspace(yb, yt, n + 1)
        mymesh(xp, yp, w.T, 'x', 'y', 'temperature')

    print("best overall temp: ", max_temp)
    """"
    We found the best overall temperature reading to be 52.945 when m = n = 100.
    We feel that this is because there is less error when the partition is smaller,
    however we have to cap out our partition somewhere so we chose to do that at 100.
    """
    return


def sa3():
    """
    Find the maximum power that can be dissipated by a 4×4cm fin while keeping the
    maximum temperature less than 80◦C. Assume the power input is the same as in
    steps 1 and 2.
    """
    xl = 0
    xr = 4
    yb = 0
    yt = 4
    u_b = 20

    power = 5.666  # watts
    length = 4  # cm
    width = .1  # cm
    K = 1.68
    H = 0.005
    m = 100
    n = 100
    w = poisson(xl, xr, yb, yt, m, n, power, length, width, K, H) + u_b

    while max(w[0]) <= 80.0:
        power += .0001
        w = poisson(xl, xr, yb, yt, m, n, power, length, width, K, H) + u_b
        # if (max(w[0]) > 80):
        #     power -= 0.05
        #     w = poisson(xl, xr, yb, yt, m, n, power, length, width, K, H)
        print('power: ', power)
        print('temp: ', max(w[0]))
    # np.set_print_options(precision=8,line_width=140)
    xp = np.linspace(xl, xr, m + 1)
    yp = np.linspace(yb, yt, n + 1)
    mymesh(xp, yp, w.T, 'x', 'y', 'w')
    """
    We found the maximum power that can be dissipated by a 4x4cm fan while maintaining a maximum
    temperature less than 80◦C was 5.666 watts.  The highest temperature that our cooling fin 
    reached given these conditions and at this power was 79.997◦C.  To find this value we ran 
    an iterative process that would recalculate the maximum temperature given a power, and would 
    increment the power value and recalculate the maximum temperature reached.  This iterative 
    process terminated once the temperature exceeded 80◦C.
    """
    return


def sa4():
    """
    Replace the aluminum fin with a copper fin, with thermal conductivity
    K=3.85W/cm◦C. Find the maximum power that can be dissipated by a 4×4 cm fin with
    the 2cm power input as in steps 1 and 2, while keeping the maximum temperature
    below 80◦C.
    """
    xl = 0
    xr = 4
    yb = 0
    yt = 4
    u_b = 20

    power = 7.149  # watts
    length = 4  # cm
    width = .1  # cm
    K = 3.85
    H = 0.005
    m = 100
    n = 100
    w = poisson(xl, xr, yb, yt, m, n, power, length, width, K, H) + u_b

    while max(w[0]) <= 80.0:
        power += .001
        w = poisson(xl, xr, yb, yt, m, n, power, length, width, K, H) + u_b
        print('power: ', power)
        print('temp: ', max(w[0]))
    power -= .001
    xp = np.linspace(xl, xr, m + 1)
    yp = np.linspace(yb, yt, n + 1)
    mymesh(xp, yp, w.T, 'x', 'y', 'w')
    """
    After replacing the aluminum fin with a copper fin, we found the maximum power that can 
    be dissipated by a 4x4cm fan while maintaining a maximum temperature less than 80◦C was 
    7.150 watts.  The highest temperature that our cooling fin reached given these conditions 
    and at this power was 79.997◦C.  To find this value we ran the same iterative process 
    as in SA3.
    """
    return


def sa5():
    """
    Plot the maximum power that can be dissipated in step 4 (keeping maximum
    temperature below 80 degrees) as a function of thermal conductivity,
    for 1.0 ≤ K ≤ 5 W/cm◦C.
    """
    xl = 0
    xr = 4
    yb = 0
    yt = 4
    u_b = 20

    power = 0.0  # watts
    length = 4  # cm
    width = .1  # cm
    x = np.linspace(1.0, 5.0, 50)
    y = np.array([])
    H = 0.005
    m = 40
    n = 40

    for i in x:
        K = i
        w = poisson(xl, xr, yb, yt, m, n, power, length, width, K, H) + u_b
        while max(w[0]) <= 80.0:
            power += .1
            w = poisson(xl, xr, yb, yt, m, n, power, length, width, K, H) + u_b
            print('power: ', power)
            print('temp: ', max(w[0]))
        power -= 0.1
        y = np.append(y, power)

    plt.plot(x, y, 'r-', linewidth=1.0, label='test')
    plt.xlabel('K (Thermal Conductivity) (W/cm◦C)')
    plt.ylabel('Max Power (Watts)')
    plt.title('SA 5')
    plt.show()

    return


def sa6():
    """
    Redo step 4 for a water-cooled fin. Assume that water has a convective heat
    transfer coefficient of H=0.1Wcm2◦C, and that the ambient water temperature
    is maintained at 20◦C
    """
    xl = 0
    xr = 4
    yb = 0
    yt = 4
    u_b = 20.0

    power = 36.515  # watts
    length = 4     # cm
    width = .1     # cm
    K = 3.85
    H = 0.1
    m = 100
    n = 100
    w = poisson(xl, xr, yb, yt, m, n, power, length, width, K, H) + u_b

    while max(w[0]) <= 80.0:
        power += 0.001
        w = poisson(xl, xr, yb, yt, m, n, power, length, width, K, H) + u_b
        print('power: ', power)
        print('temp: ', max(w[0]))
    # Since we went above 80.0, we must go back to the closest to 80.0 that we accomplished
    power -= 0.001
    xp = np.linspace(xl, xr, m + 1)
    yp = np.linspace(yb, yt, n + 1)
    mymesh(xp, yp, w.T, 'x', 'y', 'w')
    """
    Under a water cooled fin, we found the maximum power that can be dissipated by a 4x4cm fan 
    while maintaining a maximum temperature less than 80◦C was 36.517 watts.  The highest 
    temperature that our cooling fin reached given these conditions and at this power was 80.000◦C.  
    To find this value we ran the same iterative process as in SA3.
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
    xr = 4
    yb = 0
    yt = 4
    u_b = 20.0

    power = 5  # watts
    length = 4  # cm
    width = .1  # cm
    K = 1.68
    H = 0.005
    m = 10
    n = 10

    # default power source location from SA2
    w = poisson(xl, xr, yb, yt, m, n, power, length, width, K, H) + u_b
    # finds the max of a 2 dimensional array
    max_w = max(w.max(axis=1))
    print('max w: ', max_w)
    xp = np.linspace(xl, xr, m + 1)
    yp = np.linspace(yb, yt, n + 1)
    mymesh(xp, yp, w.T, 'x (cm)', 'y (cm)', 'temperature (◦C)', 'default power source location [0.0, 2.0]')

    # close to optimal power source location
    p_int = [1.25, 3.25]
    shifted_w = poisson(xl, xr, yb, yt, m, n, power, length, width, K, H, p_int) + u_b
    # finds the max of a 2 dimensional array
    max_shifted_w = max(shifted_w.max(axis=1))
    print('other max shifted w: ', max_shifted_w)
    xp = np.linspace(xl, xr, m + 1)
    yp = np.linspace(yb, yt, n + 1)
    mymesh(xp, yp, shifted_w.T, 'x (cm)', 'y (cm)', 'temperature (◦C)', 'close to optimal power source location [1.25, 3.25]')

    # optimal power source location
    p_int = [1.0, 3.0]
    shifted_w = poisson(xl, xr, yb, yt, m, n, power, length, width, K, H, p_int) + u_b
    # finds the max of a 2 dimensional array
    max_shifted_w = max(shifted_w.max(axis=1))
    print('max shifted w: ', max_shifted_w)
    xp = np.linspace(xl, xr, m + 1)
    yp = np.linspace(yb, yt, n + 1)
    mymesh(xp, yp, shifted_w.T, 'x (cm)', 'y (cm)', 'temperature (◦C)', 'optimal power source location [1.0, 3.0]')

    """
    After much trial and error, we found that the optimal location for the power
    source to be applied is right in the middle of the left side.  We found that 
    this location optimizes the dissipation across the entire cooling fin, thus
    minimizing the maximum temperature found on the cooling fin.
    """

    return

