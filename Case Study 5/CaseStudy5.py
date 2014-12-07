import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from time import clock


def K(z):
    return 10E-8 * np.exp(z / 5.)


def gamma(z):
    return 1. - ((z - 40.) / 10.) ** 2 + (1. / 2.) * ((z - 40.) / 10.) ** 4


def R(y_1, y_2, t):
    '''
    Find the reaction rate of the system at state c and time t.
    '''

    y_3 = 3.7e16          # Concentration of O_2 (constant)

    k_1 = 1.63E-16       # Reaction rate [O + O_2 ->     O_3]
    k_2 = 4.66E-16       # Reaction rate [O + O_3 -> 2 * O_2]

    a_3 = 22.62          # Constant used in calculation of k_3
    a_4 = 7.601          # Constant used in calculation of k_4
    w = np.pi / 43200    # Cycle (half a day) [1/sec]

    if np.sin(w * t) > 0:
        k_3 = np.exp(-a_3 / np.sin(w * t))
        k_4 = np.exp(-a_4 / np.sin(w * t))
    else:
        k_3 = 0
        k_4 = 0

    R_1 = -k_1 * y_1 * y_3 - k_2 * y_1 * y_2 + 2. * k_3 * y_3 + k_4 * y_2
    R_2 = +k_1 * y_1 * y_3 - k_2 * y_1 * y_2 - k_4 * y_2
    return R_1, R_2


def system(t, y):
    f = np.zeros(len(y))

    R1, R2 = R(y[0], y[1], t)
    f[0] = (dz ** -2 * (K(3. / 2.) * y[2] - (K(3. / 2.) + K(1. / 2.)) * y[0] + K(1. / 2.) * y[2]) + R1)
    f[1] = (dz ** -2 * (K(3. / 2.) * y[3] - (K(3. / 2.) + K(1. / 2.)) * y[1] + K(1. / 2.) * y[3]) + R2)

    for i in range(1, M):
        R1, R2 = R(y[2 * i], y[2 * i + 1], t)
        f[2 * i] =     (dz ** -2 * (K(i + 1. + 1. / 2.) * y[2 * i + 2] - (K(i + 1. + 1. / 2.) + K(i + 1 - 1. / 2.)) * y[2 * i]     + K(i + 1 - 1. / 2.) * y[2 * i - 2]) + R1)
        f[2 * i + 1] = (dz ** -2 * (K(i + 1. + 1. / 2.) * y[2 * i + 3] - (K(i + 1. + 1. / 2.) + K(i + 1 - 1. / 2.)) * y[2 * i + 1] + K(i + 1 - 1. / 2.) * y[2 * i - 1]) + R2)

    R1, R2 = R(y[2 * M - 2], y[2 * M - 1], t)
    f[-2] = (dz ** -2 * (K(M - 1. / 2.) * y[2 * M - 4] - (K(M - 1. / 2.) + K(M - 3. / 2.)) * y[2 * M - 2] + K(M - 3. / 2.) * y[2 * M - 4]) + R1)
    f[-1] = (dz ** -2 * (K(M - 1. / 2.) * y[2 * M - 3] - (K(M - 1. / 2.) + K(M - 3. / 2.)) * y[2 * M - 1] + K(M - 3. / 2.) * y[2 * M - 3]) + R2)

    return f


def jacobian(t, y):
    main  = np.zeros(len(y))
    sub_1 = np.zeros(len(y))
    sub_2 = np.zeros(len(y))
    sup_1 = np.zeros(len(y))
    sup_2 = np.zeros(len(y))

    y_3 = 3.7e16          # Concentration of O_2 (constant)

    k_1 = 1.63E-16       # Reaction rate [O + O_2 ->     O_3]
    k_2 = 4.66E-16       # Reaction rate [O + O_3 -> 2 * O_2]

    a_4 = 7.601          # Constant used in calculation of k_4
    w = np.pi / 43200    # Cycle (half a day) [1/sec]

    if np.sin(w * t) > 0:
        k_4 = np.exp(-a_4 / np.sin(w * t))
    else:
        k_4 = 0

    for i in range(0, M + 1, 2):
        main[2 * i] =     -dz ** -2 * (K(i + 1. + 1./2.) + K(i + 1 - 1./2.)) - k_1 * y_3 - k_2 * y[2 * i + 1]
        main[2 * i + 1] =                                                                - k_2 * y[2 * i]     - k_4

    for i in range(1, M + 1):
        sup_1[2 * i]     = -k_2 * y[2 * i] + k_4
        sup_1[2 * i + 1] = 0

    for i in range(0, M + 1):
        sub_1[2 * i]     = k_1 * y_3 - k_2 * y[2 * i + 1]
        sub_1[2 * i + 1] = 0

    sup_2[-2] = dz ** -2 * (K(i + 1. + 1./2.) + K(i + 1 - 1./2.))
    sup_2[-1] = dz ** -2 * (K(i + 1. + 1./2.) + K(i + 1 - 1./2.))
    for i in range(2, 2 * M):
        sup_2[i] = dz ** -2 *  K(i + 1. + 1./2.)

    sub_2[0:2] = dz ** -2 *  K(i + 1. - 1./2.)
    for i in range(2, 2 * M):
        sub_2[i] = dz ** -2 * (K(i + 1. + 1./2.) + K(i + 1 - 1./2.))

    jac = np.vstack((sup_2, sup_1, main, sub_1, sub_2))
    print(jac)


# Basic parameters
M = 50
dz = 20. / M
z = [30. + j * dz for j in range(M + 1)]

# This generates the initial conditions
c = np.zeros(len(2*z))
for i, _ in enumerate(z):
    c[2 * i] = 1E6 * gamma(z[i])
    c[2 * i + 1] = 1E12 * gamma(z[i])

# Time array
times = 3600. * np.array([0., 2., 4., 6., 7., 9., 12., 18., 24.])
dt = 120.

# Create result arrays
c1, c2 = [], []

# Set up ODE solver
solver = ode(system)

jacobian(0, c)


for i in range(0, len(times) - 1):
    # Initial and final time
    t_0 = times[i]
    t_f = times[i + 1]

    # Solver setup
    t = []
    sol = []
    solver.set_initial_value(c, t_0)
    # solver.set_integrator('dop853')
    # solver.set_integrator('dopri5')
    solver.set_integrator('vode', method='bdf')
    # with_jacobian=True
    start_time = clock()
    while solver.successful() and solver.t < t_f:
        solver.integrate(solver.t + dt)
        t.append(solver.t)
        sol.append(solver.y)
    print clock() - start_time, "seconds process time"

    t = np.array(t)
    sol = np.array(sol)

    # Set c1, c2 values
    y_sol = sol[-1]
    c1.append(y_sol[0::2])
    c2.append(y_sol[1::2])

    #Update IC for next iteration
    c = y_sol

fig = plt.figure()
ax1 = fig.add_subplot(2, 1, 1)
for solution in c1:
    plt.plot(z, solution)
plt.ylabel('$c_1$')
ax1.set_xticklabels([])

ax2 = fig.add_subplot(2, 1, 2)
for solution in c2:
    plt.plot(z, solution)
plt.ylabel('$c_2$')
plt.xlabel('z (km)')
plt.show()
