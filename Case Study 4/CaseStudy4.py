import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy.optimize import curve_fit
import pylab


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


def generate_things(C, s):
    constants = Config(C, s)

    # Initial Condition with boundary conditions
    T_initial = [np.sin(constants.k * x_i) for x_i in constants.x]

    # Analytic Solution
    T_analytic = Analytic(constants)
    plt.plot(T_analytic, label='Analytic')

    # Explicit Numerical Solution
    T_explicit = Explicit(T_initial, constants)
    plt.plot(T_explicit, label='Explicit')
    plt.legend()
    # plt.show()

    err = T_explicit - T_analytic
    RMS = np.sqrt(np.mean(np.square(err)))

    print('dx: ' + str(constants.dx) + ', dt: ' + str(constants.dt))
    return [constants.dx, constants.dt, RMS]


def Analytic(config):
    k, D, u, tau, x = config.k, config.D, config.u, config.tau, config.x

    N = len(x)
    t = 1 / tau
    T = np.array(x)

    for i in range(0, N):
        T[i] = np.exp(-k ** 2 * D * t) * np.sin(k * (x[i] - u * t))

    return np.array(T)


def Explicit(T, config):
    """
    FTCS (Explicit) - Forward-Time and central differencing for both the
    convective flux and the diffusive flux.
    """

    D, dt, dx, u, tau, k, L = config.D, config.dt, config.dx, config.u, config.tau, config.k, config.L

    N = len(T)
    T = np.array(T)
    T_old = np.array(T)

    t = 0
    while t < 1 / tau:
        T[0] = np.sin(k * 0)
        T[-1] = np.sin(k * L)
        for i in range(1, N - 1):
            T[i] = ((1 - 2 * D * dt / dx ** 2) * T_old[i] +
                    (D * dt / dx ** 2 + u * dt / (2 * dx)) * T_old[i - 1] +
                    (D * dt / dx ** 2 - u * dt / (2 * dx)) * T_old[i + 1])

        T_old = np.array(T)
        t += dt

    return np.array(T)


def main():
    C = [0.1,   0.5,   2, 0.5, 0.5]
    s = [0.25, 0.25, .25, 0.5,   1]

    results = []
    for C_i, s_i in zip(C, s):
        dx, dt, RMS = generate_things(C_i, s_i)
        # print('C: ' + str(C_i) + ', s: ' + str(s_i))
        results.append([dx, dt, RMS])

    print(results)
    results.sort(key=lambda x: x[0])
    print(results)
    results = np.array(results)

    plt.subplot(2, 1, 1)
    plt.plot(results[:, 0], results[:, 2], '.')
    plt.title('dx vs RMS')
    plt.yscale('log')

    logx = scipy.log10(results[:, 0])
    logy = scipy.log10(results[:, 2])

    # define our (line) fitting function
    def linear_fit(x, a, b):
        return a + b * x

    out = curve_fit(linear_fit, logx, logy)
    fit = linear_fit(results[:, 0], out[0][0], out[0][1])
    # pfinal = out[0]
    # covar = out[1]
    # print(pfinal, covar)
    plt.plot(results[:, 0], fit)
    print(out)

    plt.subplot(2, 1, 2)
    plt.plot(results[:, 1], results[:, 2], '.')
    plt.title('dt vs RMS')
    plt.yscale('log')
    plt.show()


if __name__ == "__main__":
    main()
