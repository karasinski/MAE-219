import numpy as np
import matplotlib.pyplot as plt


def Explicit():
    # Problem Parameters
    L = 1.            # Domain lenghth       [n.d.]
    T0 = 0.           # Initial temperature  [n.d.]
    T1 = 1.           # Boundary temperature [n.d.]
    t_end = 0.3
    s = 1. / 6.
    N = 21

    # Set-up Mesh
    x = np.linspace(0, L, N)
    dx = x[1] - x[0]

    # Calculate time-step
    dt = s * dx ** 2.0
    time = 0.

    # Initial Condition
    Tnew = [T0] * N

    # Boundary conditions
    Tnew[0] = T1
    Tnew[N - 1] = T1

    Told = Tnew

    # plt.axis([0, L, T0, T1])
    plt.xlabel('Length [nd]')
    plt.ylabel('Temperature [nd]')

    # Numerical Solution
    while time <= t_end:
        for i in range(1, N - 1):
            Tnew[i] = s * Told[i + 1] +  (1 - 2.0 * s) * Told[i] + s * Told[i - 1]

        time += dt
        Told = Tnew

    # Analytical Solution
    Tanalytical = Tnew
    for i in range(1, N - 1):
        Tanalytical[i] = Analytic(x[i], time)

    # Plot Numerical Solution
    plt.title('Solution at time, t = ' + str(time)[:3])
    plt.plot(x, Tnew, linewidth=1, label='Numerical Solution')
    plt.plot(x, Tanalytical, linewidth=1, label='Analytic Solution')
    plt.legend(loc='best')
    plt.show()


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
        term = (4. / ((2. * k - 1.) * np.pi)) * \
                np.sin((2. * k - 1.) * np.pi * x) * \
                np.exp(-(2. * k - 1.) ** 2. * np.pi ** 2. * t)

        # If subtracting the term from the result doesn't change the result
        # then we've hit the point at which subtracting more terms doesn't
        # matter. Else we subtract and continue.
        # print '{0} {1}, {2:.15f}'.format(k, term, result)
        if result - term == result:
            return result
        else:
            result -= term

Explicit()
