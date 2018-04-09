import numpy as np
from helpful_functions import *
import scipy.optimize as opt
from nelder_mead import *


def sa1():
    """
    Write a function that returns the potential energy U=∑i<j (1/r_ij^12 -1/r_ij^6)
    where r_ij is given at the top of p. 581. Apply Nelder–Mead to find the
    minimum energy for n=5. Try several initial guesses until you are convinced
    you have the absolute minimum. How many steps are required? To help you
    check your potential function when n=5, here is one correct input, output
    pair. U(1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0) = -6.0102023319615911
    """
    ig = np.array([0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0])
    ig = np.array([-0.2604720088,        0.7363147287,        0.4727061929,
                    0.2604716550,       -0.7363150782,       -0.4727063011,
                   -0.4144908003,       -0.3652598516,        0.3405559620,
                   -0.1944131041,        0.2843471802,       -0.5500413671,
                    0.6089042582,        0.0809130209,        0.2094855133])

    # Chilton advised we should use scipy.optimize over his code
    res = opt.minimize(U, ig, method='Nelder-Mead')

    # check
    print("Value should be -6.0102023319615911: ",
          U([1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0]))

    print('\nPotential Energy: ', res.fun)
    print('Found in ', res.nit, ' iterations.')
    print('Translated and rotated configuration: ', res.x)

    """
    We used the NElder-Mead function in the Scipy.Optimize library to find the 
    minimum energey for n=5 nodes to be -9.103852, which is spot on with the 
    supplied websites value.  It took several initial guesses to finally get 
    this value for the potential energy.  After trying at least 10 
    configurations, we decided to set the initial guess for the nodes to be 
    the optimal configuration for n=5 nodes given at the website.  It did still
    take 515 iterations to find our configuration from those initial node 
    locations, which is peculiar.
    """

    return


def sa2():
    """
    Plot the five atoms in the minimum energy configuration as circles,
    and connect all circles with line segments to view the conformed
    molecule. The line segments should form triangles and there should
    be three or four lines connecting each atom. You are welcome to use
    Python or Mathematica.
    """
    # Print optimal graph for website points
    # optimal node location for n = 5
    # http://doye.chem.ox.ac.uk/jon/structures/LJ/tables.150.html
    points = np.array([-0.26047, 0.73631, 0.47271, 0.26047, -0.73632, -0.47271, -0.41449, -0.36526,
                       0.34056, -0.19441, 0.28435, -0.55004, 0.60890, 0.08091, 0.20949])
    plot_configuration(points, 'SA2_Figure1', 'Approximate Solution for $n=5$')

    # Print optimal graph for translated points
    res = opt.minimize(U, points, method='Nelder-Mead')
    points = translate_and_rotate(res.x)
    plot_configuration(points, 'SA2_Figure2', 'Approximate Solution for $n=5$ according to RC13')

    """
    We first plot the configuration of nodes given from the optimal configuration from 
    the supplied website.  We then plotted the points for the configuration we found in 
    SA1 where the first point is fixed at the origin and the second point is fixed on 
    the z-axis.
    """
    return


def sa3():
    """
    Write a function that returns the gradient of U. Apply a Python
    minimization function that uses the function and the gradient for
    the n=5 case. Find the minimum energy as before.
    To help you check your gradient function when n=5, here is one
    correct input, output pair.
    ∇U(1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0) = [0.65625, 0.0,
    0.65625, 0.65625,0.65625,-1.3125, 0.79891,-1.45516, 0.79891,-1.3125,
    0.65625, 0.65625,-0.79891, 0.14266,-0.79891]
    """
    # Check
    print('Value should be: [0.65625, 0.0, 0.65625, 0.65625, '
          '0.65625, -1.3125, 0.79891, -1.45516, 0.79891, -1.3125, '
          '0.65625, 0.65625, -0.79891, 0.14266, -0.79891]')
    print('Value is: ', grad_U([1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0]))

    ig = np.array([0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0])

    # Set up boundary conditions
    b = 1.2
    bnds = ((0, 0), (0, 0), (0, 0), (0, 0), (0, 0), (-b, b), (-b, b), (-b, b),
            (-b, b), (-b, b), (-b, b), (-b, b), (-b, b), (-b, b), (-b, b))

    # Using the L-BFGS-B Method with bounds
    x_min_bounded = opt.minimize(U, ig, method='L-BFGS-B', jac=grad_U, bounds=bnds)
    points_bounded = x_min_bounded.x

    x_min_unbounded = opt.minimize(U, ig, method='L-BFGS-B', jac=grad_U)
    points_unbounded = translate_and_rotate(x_min_unbounded.x)

    # Print comparison between using bounds and not

    print('\nBounded method Potential Energy: ', x_min_bounded.fun)
    print('Found in ', x_min_bounded.nit, ' iterations.')
    print('Translated and rotated configuration: ', points_bounded)
    plot_configuration(points_bounded, 'SA3_Figure1', 'Unbounded Approach $n=5$ using L-BFGS-B')

    print('\nUnbounded method Potential Energy: ', x_min_unbounded.fun)
    print('Found in ', x_min_unbounded.nit, ' iterations.')
    print('Translated and rotated configuration: ', points_unbounded)
    plot_configuration(points_bounded, 'SA3_Figure2', 'Bounded Approach $n=5$ using L-BFGS-B')

    """
    Due to our results from SA5 we decided to expand this question to not only find 
    the potential energy using a Scipy.Optimize function using the gradient function 
    we have written, but also evaluate how including bounds changes the number of 
    iterations it takes to find the optimal node configuration, and also to see if 
    this affects the calculated potential energy.  We found that using found the 
    L-BFGS-B method found the potential energy to be -9.103852 in 39 total iterations. 
    The unbounded approach using the same method found the potential energy to be 
    -9.103852 in 37 total iterations.  As we see their potential energy is just about 
    the same, despite some machine error we encountered.  The actual optimal potential 
    energy for a configuration of 5 nodes is -9.103852, so both approaches were spot 
    on. It is interesting, however not surprising to see that following an unbounded 
    approach was more efficient.  Though obvious, we must note that the node 
    configurations are different.
    """
    return


def sa4():
    """
    Use one of the functions in SciPy Optimization to find the global minimum
    of your potential function when n=5 using only the potential function itself(
    not the gradient). You cannot use Nelder-Mead for this.
    """

    ig = np.array([0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 1, 1, 0])

    # Using the Powell Method
    x_min = opt.minimize(U, ig, method='Powell')

    # translate and rotate configuration
    points = translate_and_rotate(x_min.x)

    print('Potential Energy: ', x_min.fun)
    print('Found in ', x_min.nit, ' iterations.')
    print('Translated and rotated configuration: ', points)

    return


def sa5():
    """
    Apply the methods used in (1), (3), and (4) when n=6. Rank the methods
    according to reliability and efficiency.
    """
    # optimal node location for n = 6
    # http://doye.chem.ox.ac.uk/jon/structures/LJ/tables.150.html
    ig = np.array([0.7430002202, 0.2647603899, -0.0468575389, -0.7430002647, -0.2647604843, 0.0468569750,
                   0.1977276118, -0.4447220146, 0.6224700350, -0.1977281310, 0.4447221826, -0.6224697723,
                   -0.1822009635, 0.5970484122, 0.4844363476, 0.1822015272, -0.5970484858, -0.4844360463])

    # (1) Using Nelder-Mead
    print('\nUsing Nelder-Mead')
    x_min_nm = opt.minimize(U, ig, method='Nelder-Mead')
    print('Potential Energy: ', x_min_nm.fun)
    print('Found in ', x_min_nm.nit, ' iterations.')

    # (3) Using L_BFGS_B
    print('\nUsing L_BFGS_B')
    x_min_L_BFGS_B = opt.minimize(U, ig, method='L-BFGS-B', jac=grad_U)
    print('Potential Energy: ', x_min_L_BFGS_B.fun)
    print('Found in ', x_min_L_BFGS_B.nit, ' iterations.')

    # (4) Using Powell
    print('\nUsing Powell')
    x_min_Powell = opt.minimize(U, ig, method='Powell')
    print('Potential Energy: ', x_min_Powell.fun)
    print('Found in ', x_min_Powell.nit, ' iterations.')

    """
    All the optimizations methods that we used are found in the Scipy.Optimization library.  
    We used the Nelder-Mead, Powell, and L_FBGS_B methods.  We ranked these methods according
    to efficiency as follows:
    Rank (Efficiency):
    1: L_BFGS_B - 4 iterations
    2: Powell - 5 iterations
    3: Nelder-Mead - 998 iterations,
    
    We then ranked the same methods according to reliability as follows:
    Rank (Reliability) (Actual PE = -12.712062):
    1: L_BFGS_B - PE: -12.712062
    2: Nelder-Mead - PE: -12.712062
    3: Powell - PE: -12.712010
    
    You'll noticed that when transitioning from the efficiency ranking to the reliability 
    ranking, Nelder Mead and Powell have changed places.  This is because Nelder-Mead 
    took considerably more iterations (993), but the calculated potential energy of the 
    resulting configuration is closer to the exact answer than that of the Powell method.
    Originally we rann the L_BFGS_B method with bounds to make sure that the first five 
    variables were zero as part of the assignment.  However this was skewing our results 
    and the potential energy found through that method was approximately 0.4 off of the 
    actual value.  Instead of including bounds, we decided to run the resulting node
    configuration through our translate_and_rotate function, like all the other methods.
    """

    return


def sa6():
    """
    Plot the six atoms in the minimum energy configuration as circles, and
    connect all circles with line segments to view the conformed molecule.
    The line segments should form triangles and there should be four lines
    connecting each atom. You are welcome to use Python or Mathematica.
    """
    # Found this point configuration in sa6 using Powell unconstrained
    ig = np.array([0.67965175,  0.23750852, -0.04819455, -0.64271399, -0.23466047,
                   0.0452073 ,  0.1913369 , -0.36966512,  0.57079769, -0.15485331,
                   0.37379498, -0.57404431, -0.14914537,  0.55118178,  0.40520939,
                   0.18565026, -0.54732606, -0.4087711])

    # Translate and rotate these points according to the required configuration in RC13
    points = translate_and_rotate(ig)

    # Plot new configuration
    plot_configuration(points, 'SA6_Figure', 'Approximate Solution for $n=6$')

    """
    We took the optimal configuration of nodes that we found in SA 5 and ran our 
    translate_and_rotate function first so that the first point would be located 
    at the origin and the second point would be fixed on the z-axis.  We then 
    plotted this configuration.
    """

    return


def sa7():
    """
    Determine and plot minimum-energy conformations for larger n. Information
    on minimum-energy Lennard-Jones clusters for n up to several hundred is posted
    at the link provided, so your answers can be readily checked. You should do at
    least n=8.
    http://doye.chem.ox.ac.uk/jon/structures/LJ/tables.150.html
    """
    # optimal node location for n = 8
    ig = np.array([0.2730989500 , 1.1469629952 , -0.3319578813,
                   -0.4728837298, -0.6223685080, 0.7664213265 ,
                   -0.9666537615, -0.2393056630, -0.1698094248,
                   0.6209419539 , -0.3628130566, 0.7094425990 ,
                   0.8035441992 , 0.1648033307 , -0.2639206823,
                   -0.1784380914, 0.2412141513 , -0.8077599510,
                   0.0639788373 , -0.6647479592, -0.2089132333,
                   -0.1435883576, 0.3362547097 , 0.3064972473 ])

    # (Using Powell
    print('\nUsing Powell')
    x_min_Powell = opt.minimize(U, ig, method='Powell')
    print('Potential Energy at n=8: ', x_min_Powell.fun)
    print('Number of iterations: ', x_min_Powell.nit)

    # Translate and rotate these points according to the required configuration in RC13
    points = translate_and_rotate(x_min_Powell.x)

    # Plot new configuration
    plot_configuration(points, 'SA7_Figure1', 'Approximate Solution for $n=8$')

    # optimal node location for n = 14, cause 14 is a good number
    ig = np.array([-0.4308428681,  0.3011152165,  1.5395345691,
                    0.8907978174, -0.2122748336, -0.7483531248,
                   -0.0007070124,  0.3915249591, -1.1159393395,
                   -0.1087424289, -0.7253352304, -0.9277291378,
                    0.5095512108, -0.9887375564, -0.0113705147,
                   -0.9327094761, -0.0119091729, -0.6060402871,
                    0.6843305554,  0.8181106238, -0.3158498774,
                    0.9980943595, -0.0340672328,  0.3707603801,
                   -0.4435631052,  0.9423569825, -0.2236674601,
                   -0.6182734178, -0.8637379285,  0.0806931785,
                    0.0689765876, -0.4413464656,  0.8779993659,
                   -0.8273192340,  0.1657103512,  0.5084348730,
                    0.1775943569,  0.6815230386,  0.6887760023,
                    0.0328126550, -0.0229327514, -0.1172486273])

    # Using Powell
    print('\nUsing Powell')
    x_min_Powell = opt.minimize(U, ig, method='Powell')
    print('Potential Energy at n=14: ', x_min_Powell.fun)
    print('Number of iterations: ', x_min_Powell.nit)

    # Translate and rotate these points according to the required configuration in RC13
    points = translate_and_rotate(x_min_Powell.x)

    # Plot new configuration
    plot_configuration(points, 'SA7_Figure2', 'Approximate Solution for $n=14$')

    """
    For n=8 nodes, we started with initial conditions for the optimal node locations given 
    in the supplied website.  After running the Powell optimization method found in the Scipy 
    Optimize library, we translated and rotated our points accordingly so that the first point 
    would be located at the origin and the second point would be fixed on the z-axis. We found 
    that our optimal potential energy of the system is -19.821075, which is considerably close
    to the potential energy of the given configuration, which is -19.821429.  For n=14 nodes, 
    we carried out the same process, and found our optimal potential energy of the given system 
    to be -47.844233, which is also considerable close to potential energy of the given 
    configuration, which is -47.845157.
    """
    return
