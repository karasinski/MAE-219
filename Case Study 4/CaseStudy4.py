from PrettyPlots import *
import numpy as np
import matplotlib.pyplot as plt
from scipy import log10
from scipy.optimize import curve_fit


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


def generate_solutions(C, s):
    c = Config(C, s)

    # Initial Condition with boundary conditions
    T_initial = np.sin(c.k * c.x)

    # Analytic Solution
    T_analytic = Analytic(c)
    plt.plot(T_analytic, label='Analytic')

    # Explicit Numerical Solution
    T_explicit = Explicit(T_initial, c)
    plt.plot(T_explicit, label='Explicit')

    # Explicit Numerical Solution
    T_upwind = Upwind(T_initial, c)
    plt.plot(T_upwind, label='Upwind')

    plt.legend()
    plt.show()

    err = T_explicit - T_analytic
    RMS = np.sqrt(np.mean(np.square(err)))

    # print('dx: ' + str(c.dx) + ', dt: ' + str(c.dt) + ', RMS: ' + str(RMS))
    return [c.dx, c.dt, RMS]


def Analytic(c):
    k, D, u, tau, x = c.k, c.D, c.u, c.tau, c.x

    N = len(x)
    t = 1 / tau
    T = np.array(x)

    for i in range(0, N):
        T[i] = np.exp(-k ** 2 * D * t) * np.sin(k * (x[i] - u * t))

    return np.array(T)


def Explicit(T, c):
    """
    FTCS (Explicit) - Forward-Time and central differencing for both the
    convective flux and the diffusive flux.
    """

    D, dt, dx, u, tau = c.D, c.dt, c.dx, c.u, c.tau

    N = len(T)
    T = np.array(T)
    T_old = np.array(T)

    # spatial_stability = dx < (2 * D) / u
    # temporal_stability = dt < dx ** 2 / (2 * D)
    # print('Stability: spatial ' + str(spatial_stability) +
                  # ', temporal ' + str(temporal_stability))

    t = 0
    while t < 1 / tau:
        for i in range(1, N - 1):
            T[i] = ((1 - 2 * D * dt / dx ** 2) * T_old[i] +
                    (D * dt / dx ** 2 + u * dt / (2 * dx)) * T_old[i - 1] +
                    (D * dt / dx ** 2 - u * dt / (2 * dx)) * T_old[i + 1])

        # Enforce our periodic boundary condition
        T[-1] = ((1 - 2 * D * dt / dx ** 2) * T_old[-1] +
                (D * dt / dx ** 2 + u * dt / (2 * dx)) * T_old[-2] +
                (D * dt / dx ** 2 - u * dt / (2 * dx)) * T_old[1])
        T[0] = T[-1]

        T_old = np.array(T)
        t += dt

    return np.array(T_old)


def Upwind(T, c):
    '''
    Upwind-Finite Volume method: Explicit (forward Euler), with the convective
    flux treated using the basic upwind method and the diffusive flux treated
    using central differencing.
    '''

    D, dt, dx, u, tau = c.D, c.dt, c.dx, c.u, c.tau

    N = len(T)
    T = np.array(T)
    T_old = np.array(T)

    t = 0
    while t <= 1 / tau:
        for i in range(0, N - 1):
            T[i] = (D * dt / dx ** 2 * (T_old[i + 1] - 2 * T_old[i] + T_old[i - 1]) -
                    u * dt / (2 * dx) * (3 * T_old[i] - 4 * T_old[i - 1] + T_old[i - 2]) +
                    T_old[i])

        # Enforce our periodic boundary condition
        T[-1] = (D * dt / dx ** 2 * (T_old[0] - 2 * T_old[-1] + T_old[-2]) -
                 u * dt / (2 * dx) * (3 * T_old[-1] - 4 * T_old[-2] + T_old[-3]) +
                 T_old[-1])
        T[0] = T[-1]

        T_old = np.array(T)
        t += dt

    return np.array(T_old)


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
    # print('x order: ', order_accuracy_x, 't order: ', order_accuracy_t)

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
