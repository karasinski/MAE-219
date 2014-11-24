from PrettyPlots import *
import numpy as np
import matplotlib.pyplot as plt
from scipy import log10
from scipy.optimize import curve_fit
import scipy.sparse as sparse
import os


class Config(object):
    def __init__(self, C, s):
        # Import parameters
        self.C = C
        self.s = s

        # Problem constants
        self.L = 1.                  # m
        self.D = 0.005               # m^2/s
        self.u = 0.2                 # m/s
        self.k = 2 * np.pi / self.L  # m^-1
        self.tau = 1 / (self.k ** 2 * self.D)

        # Set-up Mesh and Calculate time-step
        self.dx = self.C * self.D / (self.u * self.s)
        self.dt = self.C * self.dx / self.u
        self.x = np.append(np.arange(0, self.L, self.dx), self.L)


def stability(c):
    C, s, D, u, dx, dt = c.C, c.s, c.D, c.u, c.dx, c.dt

    print('C = ', C, ' s = ', s)
    stable_x = dx < (2 * D) / u
    stable_t = dt < dx ** 2 / (2 * D)
    print('Explicit: spatial ' + str(stable_x) + ', temporal ' + str(stable_t))

    # stable_x = ???
    # stable_t = ???
    # print('Upwind: spatial ' + str(stable_x) + ', temporal ' + str(stable_t))

    # stable_x = ???
    stable_t = dt < dx ** 2 / D
    print('Trapezoidal: temporal ' + str(stable_t))

    # stable_x = ???
    # stable_t = ???
    # print('QUICK: spatial ' + str(stable_x) + ', temporal ' + str(stable_t))


def save_figure(analytic, solution, title):
    plt.plot(analytic, label='Analytic')
    plt.plot(solution, '.', label=title.split(' ')[0])

    # Calculate NRMS for this solution
    err = solution - analytic
    NRMS = np.sqrt(np.mean(np.square(err)))/(max(analytic) - min(analytic))

    plt.ylabel('$\Phi$')
    plt.xlabel('x')
    plt.title('C=' + title.split(' ')[1] +
              ' s=' + title.split(' ')[2] +
              ' NRMS={0:.3e}'.format(NRMS))
    plt.legend(loc='best')

    # Save plots
    save_name = title + '.pdf'
    try:
        os.mkdir('figures')
    except Exception:
        pass

    plt.savefig('figures/' + save_name, bbox_inches='tight')
    plt.clf()


def save_state(analytic, solutions, state):
    plt.plot(analytic, label='Analytic')
    for solution in solutions:
        plt.plot(solution[0], '.', label=solution[1])

    plt.ylabel('$\Phi$')
    plt.xlabel('x')

    title = 'C=' + state.split(' ')[0] + ' s=' + state.split(' ')[1]
    plt.title(title)
    plt.legend(loc='best')

    # Save plots
    save_name = title + '.pdf'
    try:
        os.mkdir('figures')
    except Exception:
        pass

    plt.savefig('figures/' + save_name, bbox_inches='tight')
    plt.clf()


def generate_solutions(C, s):
    c = Config(C, s)

    # Spit out some stability information
    stability(c)

    # Initial Condition with boundary conditions
    Phi_initial = np.sin(c.k * c.x)

    # Analytic Solution
    Phi_analytic = Analytic(c)

    # Explicit Solution
    Phi_explicit = Explicit(Phi_initial, c)

    # Upwind Solution
    Phi_upwind = Upwind(Phi_initial, c)

    # Trapezoidal Solution
    Phi_trapezoidal = Trapezoidal(Phi_initial, c)

    # QUICK Solution
    Phi_quick = QUICK(Phi_initial, c)

    # Save individual comparisons
    # save_figure(Phi_analytic, Phi_explicit, 'Explicit ' + str(C) + ' ' + str(s))
    # save_figure(Phi_analytic, Phi_upwind, 'Upwind ' + str(C) + ' ' + str(s))
    # save_figure(Phi_analytic, Phi_trapezoidal, 'Trapezoidal ' + str(C) + ' ' + str(s))
    # save_figure(Phi_analytic, Phi_quick, 'QUICK ' + str(C) + ' ' + str(s))

    # Save group comparison
    solutions = [(Phi_explicit, 'Explicit'),
                 (Phi_upwind, 'Upwind'),
                 (Phi_trapezoidal, 'Trapezoidal'),
                 (Phi_quick, 'QUICK')]
    save_state(Phi_analytic, solutions, str(C) + ' ' + str(s))

    err = Phi_explicit - Phi_analytic
    RMS = np.sqrt(np.mean(np.square(err)))

    # print('dx: ' + str(c.dx) + ', dt: ' + str(c.dt) + ', RMS: ' + str(RMS))
    return [c.dx, c.dt, RMS]


def Analytic(c):
    k, D, u, tau, x, dt = c.k, c.D, c.u, c.tau, c.x, c.dt

    N = len(x)
    Phi = np.array(x)

    t = 0
    while t < tau:
        t += dt

    for i in range(0, N):
        Phi[i] = np.exp(-k ** 2 * D * t) * np.sin(k * (x[i] - u * t))

    return np.array(Phi)


def Explicit(Phi, c):
    """
    FTCS (Explicit) - Forward-Time and central differencing for both the
    convective flux and the diffusive flux.
    """

    D, dt, dx, u, tau = c.D, c.dt, c.dx, c.u, c.tau

    N = len(Phi)
    Phi = np.array(Phi)
    Phi_old = np.array(Phi)

    t = 0
    while t < tau:
        for i in range(1, N - 1):
            Phi[i] = ((1 - 2 * D * dt / dx ** 2) * Phi_old[i] +
                      (D * dt / dx ** 2 + u * dt / (2 * dx)) * Phi_old[i - 1] +
                      (D * dt / dx ** 2 - u * dt / (2 * dx)) * Phi_old[i + 1])

        # Enforce our periodic boundary condition
        Phi[-1] = ((1 - 2 * D * dt / dx ** 2) * Phi_old[-1] +
                   (D * dt / dx ** 2 + u * dt / (2 * dx)) * Phi_old[-2] +
                   (D * dt / dx ** 2 - u * dt / (2 * dx)) * Phi_old[1])
        Phi[0] = Phi[-1]

        Phi_old = np.array(Phi)
        t += dt

    return np.array(Phi_old)


def Upwind(Phi, c):
    '''
    Upwind-Finite Volume method: Explicit (forward Euler), with the convective
    flux treated using the basic upwind method and the diffusive flux treated
    using central differencing.
    '''

    D, dt, dx, u, tau = c.D, c.dt, c.dx, c.u, c.tau

    N = len(Phi)
    Phi = np.array(Phi)
    Phi_old = np.array(Phi)

    t = 0
    while t <= tau:
        # Enforce our periodic boundary condition
        Phi[0] = (D * dt / dx ** 2 * (Phi_old[1] - 2 * Phi_old[0] + Phi_old[-1]) -
                  u * dt / (2 * dx) * (3 * Phi_old[0] - 4 * Phi_old[-1] + Phi_old[-2]) +
                  Phi_old[0])

        Phi[1] = (D * dt / dx ** 2 * (Phi_old[2] - 2 * Phi_old[1] + Phi_old[0]) -
                  u * dt / (2 * dx) * (3 * Phi_old[1] - 4 * Phi_old[0] + Phi_old[-1]) +
                  Phi_old[1])

        for i in range(2, N - 1):
            Phi[i] = (D * dt / dx ** 2 * (Phi_old[i + 1] - 2 * Phi_old[i] + Phi_old[i - 1]) -
                      u * dt / (2 * dx) * (3 * Phi_old[i] - 4 * Phi_old[i - 1] + Phi_old[i - 2]) +
                      Phi_old[i])

        Phi[-1] = (D * dt / dx ** 2 * (Phi_old[0] - 2 * Phi_old[-1] + Phi_old[-2]) -
                   u * dt / (2 * dx) * (3 * Phi_old[-1] - 4 * Phi_old[-2] + Phi_old[-3]) +
                   Phi_old[-1])

        Phi_old = np.array(Phi)
        t += dt

    return np.array(Phi_old)


def Trapezoidal(Phi, c):
    D, dt, dx, u, tau = c.D, c.dt, c.dx, c.u, c.tau

    N = len(Phi)
    Phi = np.array(Phi)
    Phi_old = np.array(Phi)

    # Create Coefficient Matrix
    upper = [-(dt * D) / (2 * dx ** 2) + dt * u / (4 * dx) for _ in range(0, N)]
    main = [1 + (dt * D / (dx ** 2)) for _ in range(0, N)]
    lower = [-(dt * D) / (2 * dx ** 2) - dt * u / (4 * dx) for _ in range(0, N)]

    data = lower, main, upper
    diags = np.array([-1, 0, 1])
    matrix = sparse.spdiags(data, diags, N, N).todense()

    # Set values for cyclic boundary conditions
    matrix[0, N - 1] = -(dt * D) / (2 * dx ** 2) - dt * u / (4 * dx)
    matrix[N - 1, 0] = -(dt * D) / (2 * dx ** 2) + dt * u / (4 * dx)

    # create blank b array
    b = np.array(Phi_old)

    t = 0
    while t <= tau:
        # Enforce our periodic boundary condition
        b[0] = ((dt * D / (2 * dx ** 2)) * (Phi_old[1] - 2 * Phi_old[0] + Phi_old[N - 1]) -
                (u * dt / (4 * dx)) * (Phi_old[1] - Phi_old[N - 1]) +
                Phi_old[0])

        for i in range(1, N - 1):
            b[i] = ((dt * D / (2 * dx ** 2)) * (Phi_old[i + 1] - 2 * Phi_old[i] + Phi_old[i - 1]) -
                    (u * dt / (4 * dx)) * (Phi_old[i + 1] - Phi_old[i - 1]) +
                    Phi_old[i])

        # Enforce our periodic boundary condition
        b[N - 1] = ((dt * D / (2 * dx ** 2)) * (Phi_old[0] - 2 * Phi_old[N - 1] + Phi_old[N - 2]) -
                    (u * dt / (4 * dx)) * (Phi_old[0] - Phi_old[N - 2]) +
                     Phi_old[N - 1])

        # Solve matrix
        Phi = np.linalg.solve(matrix, b)

        Phi_old = np.array(Phi)
        t += dt

    return np.array(Phi_old)


def QUICK(Phi, c):
    D, dt, dx, u, tau = c.D, c.dt, c.dx, c.u, c.tau

    N = len(Phi)
    Phi = np.array(Phi)
    Phi_old = np.array(Phi)

    t = 0
    while t <= tau:
        # Enforce our periodic boundary condition
        Phi[0] = (dt * D / dx ** 2 * (Phi_old[1] - 2 * Phi_old[0] + Phi_old[N - 1]) -
                  dt * u / (8 * dx) * (3 * Phi_old[1] - Phi_old[N - 2] + 6 * Phi_old[N - 1] - 8 * Phi_old[0]) +
                  Phi_old[0])
        Phi[1] = (dt * D / dx ** 2 * (Phi_old[2] - 2 * Phi_old[1] + Phi_old[0]) -
                  dt * u / (8 * dx) * (3 * Phi_old[2] - Phi_old[N - 1] + 6 * Phi_old[0] - 8 * Phi_old[1]) +
                  Phi_old[1])

        for i in range(2, N - 1):
            Phi[i] = (dt * D / dx ** 2 * (Phi_old[i + 1] - 2 * Phi_old[i] + Phi_old[i - 1]) -
                      dt * u / (8 * dx) * (3 * Phi_old[i + 1] - Phi_old[i - 2] + 6 * Phi_old[i - 1] - 8 * Phi_old[i]) +
                      Phi_old[i])

        # Enforce our periodic boundary condition
        Phi[N - 1] = (dt * D / dx ** 2 * (Phi_old[0] - 2 * Phi_old[N - 1] + Phi_old[N - 2]) -
                      dt * u / (8 * dx) * (3 * Phi_old[0] - Phi_old[N - 3] + 6 * Phi_old[N - 2] - 8 * Phi_old[N - 1]) +
                      Phi_old[N - 1])

        # Increment
        Phi_old = np.array(Phi)
        t += dt

    return np.array(Phi_old)


def linear_fit(x, a, b):
    '''Define our (line) fitting function'''
    return a + b * x


def effective_order(x, y):
    '''Find slope of log log plot to find our effective order of accuracy'''

    logx = log10(x)
    logy = log10(y)
    out = curve_fit(linear_fit, logx, logy)

    return out[0][1]


def main():
    C = [0.1,   0.5,   2, 0.5, 0.5]
    s = [0.25, 0.25, .25, 0.5,   1]

    results = []
    for C_i, s_i in zip(C, s):
        dx, dt, RMS = generate_solutions(C_i, s_i)
        # print('C: ' + str(C_i) + ', s: ' + str(s_i))
        results.append([dx, dt, RMS])

    # Sort and convert
    results.sort(key=lambda x: x[0])
    results = np.array(results)

    # Pull out data
    x = results[:, 0]
    t = results[:, 1]
    RMS = results[:, 2]

    # Find effective order of accuracy
    order_accuracy_x = effective_order(x, RMS)
    order_accuracy_t = effective_order(t, RMS)
    print('x order: ', order_accuracy_x, 't order: ', order_accuracy_t)

    # Show effect of dx on RMS
    plt.subplot(2, 1, 1)
    plt.plot(x, RMS, '.')
    plt.title('dx vs RMS, effective order {0:1.2f}'.format(order_accuracy_x))
    plt.xscale('log')
    plt.yscale('log')

    # Show effect of dt on RMS
    plt.subplot(2, 1, 2)
    plt.plot(t, RMS, '.')
    plt.title('dt vs RMS, effective order {0:1.2f}'.format(order_accuracy_t))
    plt.xscale('log')
    plt.yscale('log')

    # Finally show it off
    plt.show()


if __name__ == "__main__":
    main()
