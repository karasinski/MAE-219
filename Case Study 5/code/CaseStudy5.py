import numpy as np
from scipy.integrate import ode
from time import clock
from PrettyPlots import *


def K(z):
    zz = 30. + z * dz
    return 1E-8 * np.exp(zz / 5.)


def gamma(z):
    return 1. - ((z - 40.) / 10.) ** 2 + (1. / 2.) * ((z - 40.) / 10.) ** 4


def R(y_1, y_2, t):
    '''
    Find the reaction rates, R_1 and R_2, of the system at state c and time t.
    '''

    if np.sin(w * t) > 0.:
        k_3 = np.exp(-a_3 / np.sin(w * t))
        k_4 = np.exp(-a_4 / np.sin(w * t))
    else:
        k_3 = 0.
        k_4 = 0.

    R_1 = -k_1 * y_1 * y_3 - k_2 * y_1 * y_2 + 2. * k_3 * y_3 + k_4 * y_2
    R_2 = +k_1 * y_1 * y_3 - k_2 * y_1 * y_2 - k_4 * y_2
    return R_1, R_2


def system(t, y):
    f = np.zeros(len(y))

    R1, R2 = R(y[0], y[1], t)
    l_p, l_m = 3. / 2., 1. / 2.

    f[0] = (dz ** -2 * (K(l_p) * y[2] - (K(l_p) + K(l_m)) * y[0] + K(l_m) * y[2]) + R1)
    f[1] = (dz ** -2 * (K(l_p) * y[3] - (K(l_p) + K(l_m)) * y[1] + K(l_m) * y[3]) + R2)

    for i in range(1, M):
        R1, R2 = R(y[2 * i], y[2 * i + 1], t)
        l_p, l_m = 1. * i + 3. / 2., 1. * i + 1. / 2.

        f[2 * i] =     (dz ** -2 * (K(l_p) * y[2 * i + 2] - (K(l_p) + K(l_m)) * y[2 * i]     + K(l_m) * y[2 * i - 2]) + R1)
        f[2 * i + 1] = (dz ** -2 * (K(l_p) * y[2 * i + 3] - (K(l_p) + K(l_m)) * y[2 * i + 1] + K(l_m) * y[2 * i - 1]) + R2)

    R1, R2 = R(y[2 * M], y[2 * M + 1], t)
    l_p, l_m = 1. * M + 1. / 2., 1. * M - 1. / 2.

    f[2 * M]     = (dz ** -2 * (K(l_p) * y[2 * M - 2] - (K(l_p) + K(l_m)) * y[2 * M]     + K(l_m) * y[2 * M - 2]) + R1)
    f[2 * M + 1] = (dz ** -2 * (K(l_p) * y[2 * M - 1] - (K(l_p) + K(l_m)) * y[2 * M + 1] + K(l_m) * y[2 * M - 1]) + R2)

    return f


def solve(solver, c, time, integrator):
    # Create result arrays
    c1, c2, c1_40km, c2_40km, t = [], [], [], [], []

    start_time = clock()
    for i in range(0, len(time) - 1):
        # Initial and final time
        t_0 = time[i]
        t_f = time[i + 1]

        # Solver setup
        sol = []
        solver.set_initial_value(c, t_0)
        while solver.successful() and solver.t < t_f:
            solver.integrate(solver.t + dt)
            sol.append(solver.y)

            # keep time history for 40km point
            one, two = sol[-1][0::2], sol[-1][1::2]
            mid_one, mid_two = one[M / 2], two[M / 2]
            c1_40km.append(mid_one), c2_40km.append(mid_two)
            t.append(solver.t)

            print "{0:03.2f}%".format(100. * solver.t / time[-1])

        # Save c1, c2 solutions
        c1.append(one), c2.append(two)

        #Update initial conditions for next iteration
        c = sol[-1]

    elapsed_time = clock() - start_time
    print(elapsed_time, "seconds process time")

    output = [c1, c2, c1_40km, c2_40km, t]
    return output


def save_variables(name, z, c1, c2, t, c1_40km, c2_40km):
    try:
        os.mkdir('data')
    except Exception:
        pass

    try:
        os.mkdir('data/' + name)
    except Exception:
        pass

    np.savetxt('data/' + name + '/z.csv', z)
    np.savetxt('data/' + name + '/c1.csv', c1)
    np.savetxt('data/' + name + '/c2.csv', c2)
    np.savetxt('data/' + name + '/t.csv', t)
    np.savetxt('data/' + name + '/c1_40km.csv', c1_40km)
    np.savetxt('data/' + name + '/c2_40km.csv', c2_40km)


def load_variables(name):
    z       = np.loadtxt('data/' + name + '/z.csv')
    c1      = np.loadtxt('data/' + name + '/c1.csv')
    c2      = np.loadtxt('data/' + name + '/c2.csv')
    t       = np.loadtxt('data/' + name + '/t.csv')
    c1_40km = np.loadtxt('data/' + name + '/c1_40km.csv')
    c2_40km = np.loadtxt('data/' + name + '/c2_40km.csv')

    return z, c1, c2, t, c1_40km, c2_40km


def run_trials(z, integrators, times, M):
    # Set up ODE solver
    for integrator in integrators:
        if integrator == 'dop853' or integrator == 'dopri5':
            solver = ode(system)
            solver.set_integrator(integrator, atol=1E-1, rtol=1E-3)
        elif integrator == 'bdf':
            solver = ode(system)
            solver.set_integrator('vode', method=integrator, atol=1E-1, rtol=1E-3, nsteps=2000)

        name = integrator + ' ' + str(times[-1]) + ' ' + str(M)
        try:
            z, c1, c2, t, c1_40km, c2_40km = load_variables(name)
            print "Loaded data for: " + name
        except:
            print "Starting solver: ", integrator, "with times", times
            c1, c2, c1_40km, c2_40km, t = solve(solver, c, times, integrator)
            save_variables(name, z, c1, c2, t, c1_40km, c2_40km)

        # And plot some things
        # labels = [str(int(time / 3600.)) + " hours" for time in times[1:]]
        # plot_c1(z, c, c1, labels, name)
        # plot_c2(z, c, c2, labels, name)
        # plot_40km(t, c1_40km, c2_40km, name)


def sensitivity_analysis(integrators, times, meshes):
    plt.figure()
    for integrator in integrators:

        z_M, c1_M, c2_M = [], [], []
        for M in meshes:
            name = integrator + ' ' + str(times[-1]) + ' ' + str(M)
            try:
                z, c1, c2, _, _, _ = load_variables(name)
            except Exception:
                print Exception
            z_M.append(list(z))
            c1_M.append(list(c1[-1]))
            c2_M.append(list(c2[-1]))

        best_z = z_M[-1]
        best_c1, best_c2 = c1_M[-1], c2_M[-1]
        NRMS1, NRMS2 = [], []
        for j, mesh in enumerate(z_M):
            if j + 1 == len(z_M): break     # RMS with yourself is silly
            best1, best2, curr1, curr2 = [], [], [], []
            for i, element in enumerate(best_z):
                if element in mesh:
                    best1.append(best_c1[i])
                    best2.append(best_c2[i])
                    curr1.append(c1_M[j][mesh.index(element)])
                    curr2.append(c2_M[j][mesh.index(element)])

            best1, best2 = np.array(best1), np.array(best2)
            curr1, curr2 = np.array(curr1), np.array(curr2)

            err1, err2 = curr1 - best1, curr2 - best2
            NRMS1.append(np.sqrt(np.mean(np.square(err1)))/(max(best1) - min(best1)))
            NRMS2.append(np.sqrt(np.mean(np.square(err2)))/(max(best2) - min(best2)))
            # print meshes[j], NRMS1, NRMS2

        plt.plot([mesh for mesh in meshes][0:-1], NRMS1, ':', label= integrator + ' $c_1$')
        plt.plot([mesh for mesh in meshes][0:-1], NRMS2, '--', label= integrator + ' $c_2$')

    plt.ylabel('NRMS')
    plt.xlabel('M')
    plt.xlim([meshes[0] - 1, meshes[-2] + 1])
    plt.legend()
    save_name = str(meshes) + '.pdf'
    save_plot(save_name)


# Basic problem parameters
y_3 = 3.7E16         # Concentration of O_2 (constant)
k_1 = 1.63E-16       # Reaction rate [O + O_2 ->     O_3]
k_2 = 4.66E-16       # Reaction rate [O + O_3 -> 2 * O_2]
a_3 = 22.62          # Constant used in calculation of k_3
a_4 = 7.601          # Constant used in calculation of k_4
w = np.pi / 43200.   # Cycle (half a day) [1/sec]
dt = 60.

# Base Case
M = 50               # Number of sections
dz = 20. / M         # 20km divided by M subsections

# This generates the initial conditions
c = np.zeros(2 * (M + 1))
z = np.zeros(M + 1)
for i in range(0, M + 1):
    z[i] = 30. + i * dz
    c[2 * i]     = 1E6  * gamma(z[i])
    c[2 * i + 1] = 1E12 * gamma(z[i])

# Run the trials
integrators = ['dopri5']
times = 3600. * np.array([0., 2., 4., 6., 7., 9., 12., 18., 24.])
run_trials(z, integrators, times)

integrators = ['dopri5', 'bdf']
times = 3600. * np.array([0., 2., 4., 6., 7., 9., 12., 18., 240.])
run_trials(z, integrators, times)

# Mesh Analysis
meshes = [5, 10, 20, 40, 50, 80, 160]
for M in meshes:
    dz = 20. / M         # 20km divided by M subsections

    # This generates the initial conditions
    c = np.zeros(2 * (M + 1))
    z = np.zeros(M + 1)
    for i in range(0, M + 1):
        z[i] = 30. + i * dz
        c[2 * i]     = 1E6  * gamma(z[i])
        c[2 * i + 1] = 1E12 * gamma(z[i])

    # Time array
    dt = 60.

    integrators = ['dopri5', 'bdf']
    times = 3600. * np.array([0., 2., 4.])
    run_trials(z, integrators, times, M)

sensitivity_analysis(integrators, times, meshes)
