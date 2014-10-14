import numpy as np
import matplotlib.pyplot as plt


def Solver(s, time, show_plot=False):
    # Problem Parameters
    L = 1.            # Domain lenghth       [n.d.]
    T0 = 0.           # Initial temperature  [n.d.]
    T1 = 1.           # Boundary temperature [n.d.]
    N = 21

    # Set-up Mesh
    x = np.linspace(0, L, N)
    dx = x[1] - x[0]

    # Calculate time-step
    dt = s * dx ** 2.0

    # Initial Condition with boundary conditions
    T_initial = [T0] * N
    T_initial[0] = T1
    T_initial[N - 1] = T1

    # Explicit Numerical Solution
    T_explicit = Explicit(T_initial, time, dt, s)

    # Implicit Numerical Solution
    T_implicit = Implicit(T_initial, time, dt, s)

    # Analytical Solution
    T_analytic = [0] * N
    for i in range(0, N):
        T_analytic[i] = Analytic(x[i], time)

    # Find the RMS
    RMS = RootMeanSquare(T_explicit, T_analytic)

    # plt.axis([0, L, T0, T1])
    plt.xlabel('Length [nd]')
    plt.ylabel('Temperature [nd]')

    # Plot Numerical Solution
    num = '{:.2e}'.format(RMS)
    plt.title('s = ' + str(s)[:5] + ', t = ' + str(time)[:3] + ', RMS = ' + num)
    plt.plot(x, T_explicit, 'xr', linewidth=1, label='Explicit Solution')
    plt.plot(x, T_implicit, '+g', linewidth=1, label='Implicit Solution')
    plt.plot(x, T_analytic, 'ob', mfc='none', linewidth=1, label='Analytic Solution')
    plt.legend(loc='best')

    save_name = 'proj_1_s_' + str(s)[:5] + '_t_' + str(time) + '.png'
    plt.savefig(save_name, bbox_inches='tight')
    if show_plot:
        plt.show()
    plt.clf()

    return RMS


def Explicit(Told, t_end, dt, s):
    N = len(Told)
    time = 0.
    Tnew = Told

    while time <= t_end:
        for i in range(1, N - 1):
            Tnew[i] = s * Told[i + 1] + (1 - 2.0 * s) * Told[i] + s * Told[i - 1]

        time += dt
        Told = Tnew

    return Told


def Implicit(Told, t_end, dt, s):
    N = len(Told)
    time = 0.

    # Build our 'A' matrix
    a = [-s] * N
    a[0] = 0
    a[-1] = 0
    b = [1 + 2*s] * N
    b[0] = 1            # hold boundary
    b[-1] = 1           # hold boundary
    c = a

    while time <= t_end:
        Tnew = TDMAsolver(a, b, c, Told)

        time += dt
        Told = Tnew

    return Told


def RootMeanSquare(a, b):
    N = len(a)

    RMS = 0.
    for i in range(0, N):
        RMS += (a[i] - b[i])**2.

    RMS = RMS ** (1./2.)
    RMS /= N

    return RMS


def TDMAsolver(a, b, c, d):
    '''
    Tridiagonal Matrix Algorithm (a.k.a Thomas algorithm).

    a b c d can be NumPy array type or Python list type.
    '''

    nf = len(a)     # number of equations
    ac, bc, cc, dc = map(np.array, (a, b, c, d))     # copy the array

    for it in range(1, nf):
        mc = ac[it] / bc[it - 1]
        bc[it] = bc[it] - mc * cc[it - 1]
        dc[it] = dc[it] - mc * dc[it - 1]

    xc = ac
    xc[-1] = dc[-1] / bc[-1]

    for il in range(nf - 2, -1, -1):
        xc[il] = (dc[il] - cc[il] * xc[il + 1]) / bc[il]

    del bc, cc, dc  # delete variables from memory

    return xc


def Analytic(x, t):
    # The analytic answer is 1 - Sum(terms). Though there are an infinite
    # number of terms, only the first few matter when we compute the answer.

    result = 1
    large_number = 1E6

    for k in range(1, int(large_number) + 1):
        term = ((4. / ((2. * k - 1.) * np.pi)) *
                np.sin((2. * k - 1.) * np.pi * x) *
                np.exp(-(2. * k - 1.) ** 2. * np.pi ** 2. * t))

        # If subtracting the term from the result doesn't change the result
        # then we've hit the point at which subtracting more terms doesn't
        # matter. Else we subtract and continue.
        # print '{0} {1}, {2:.15f}'.format(k, term, result)
        if result - term == result:
            return result
        else:
            result -= term


def main():
    # Loop over requested values for s and t
    s = [1. / 6., .25, .5, .75]
    t = [0.3, 0.6, 0.9]
    RMS = []
    for i, s_ in enumerate(s):
        sRMS = [0] * len(t)
        for j, t_ in enumerate(t):
            sRMS[j] = Solver(s_, t_, False)
            # print i, j, sRMS[j]
        RMS.append(sRMS)

    plt.figure()
    plt.plot(t, RMS[0], '.r', label='s = 1/6')
    plt.plot(t, RMS[1], '.g', label='s = .25')
    plt.plot(t, RMS[2], '.b', label='s = .50')
    plt.plot(t, RMS[3], '.k', label='s = .75')
    plt.xlabel('t')
    plt.ylabel('RMS')
    plt.title('RMS vs t')
    plt.legend(loc='best')

    save_name = 'proj_1_rms_vs_t.png'
    plt.savefig(save_name, bbox_inches='tight')
    # plt.show()
    plt.clf()

    # Convert to np array to make this easier...
    RMS = np.array(RMS)

    plt.figure()
    plt.plot(s, RMS[:, 0], '.r', label='t = 0.3')
    plt.plot(s, RMS[:, 1], '.g', label='t = 0.6')
    plt.plot(s, RMS[:, 2], '.b', label='t = 0.9')
    plt.xlabel('s')
    plt.ylabel('RMS')
    plt.title('RMS vs s')
    plt.legend(loc='best')

    save_name = 'proj_1_rms_vs_s.png'
    plt.savefig(save_name, bbox_inches='tight')
    # plt.show()
    plt.clf()

if __name__ == "__main__":
    main()
