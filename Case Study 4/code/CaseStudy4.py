from PrettyPlots import *
import numpy as np
import scipy.sparse as sparse


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


def Analytic(c):
    k, D, u, tau, x = c.k, c.D, c.u, c.tau, c.x

    N = len(x)
    Phi = np.array(x)

    for i in range(0, N):
        Phi[i] = np.exp(-k ** 2 * D * tau) * np.sin(k * (x[i] - u * tau))

    return np.array(Phi)


def FTCS(Phi, c):
    """
    FTCS (Explicit) - Forward-Time and central differencing for both the
    convective flux and the diffusive flux.
    """

    D, dt, dx, u, tau = c.D, c.dt, c.dx, c.u, c.tau

    N = len(Phi)
    Phi = np.array(Phi)
    Phi_old = np.array(Phi)

    A = D * dt / dx ** 2
    B = u * dt / (2 * dx)

    t = 0
    while t < tau:
        for i in range(1, N - 1):
            Phi[i] = ((A + B) * Phi_old[i - 1] +
                      (1 - 2 * A) * Phi_old[i] +
                      (A - B) * Phi_old[i + 1])

        # Enforce our periodic boundary condition
        Phi[-1] = ((A + B) * Phi_old[-2] +
                   (1 - 2 * A) * Phi_old[-1] +
                   (A - B) * Phi_old[1])
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

    A = D * dt / dx ** 2
    B = u * dt / (2 * dx)

    t = 0
    while t <= tau:
        for i in range(2, N - 1):
            Phi[i] = (A * (Phi_old[i + 1] - 2 * Phi_old[i] + Phi_old[i - 1]) -
                      B * (3 * Phi_old[i] - 4 * Phi_old[i - 1] + Phi_old[i - 2]) +
                      Phi_old[i])

        Phi[-1] = (A * (Phi_old[1] - 2 * Phi_old[-1] + Phi_old[-2]) -
                   B * (3 * Phi_old[-1] - 4 * Phi_old[-2] + Phi_old[-3]) +
                   Phi_old[-1])
        Phi[0] = Phi[-1]
        Phi[1] = (A * (Phi_old[2] - 2 * Phi_old[1] + Phi_old[-1]) -
                  B * (3 * Phi_old[1] - 4 * Phi_old[-1] + Phi_old[-2]) +
                  Phi_old[1])

        Phi_old = np.array(Phi)
        t += dt

    return np.array(Phi_old)


def Trapezoidal(Phi, c):
    D, dt, dx, u, tau = c.D, c.dt, c.dx, c.u, c.tau

    N = len(Phi)
    Phi = np.array(Phi)
    Phi_old = np.array(Phi)

    A = dt * D / (2 * dx**2)
    B = dt * u / (4 * dx)

    # Create Coefficient Matrix
    lower = [-A - B for _ in range(0, N)]
    main = [1 + 2 * A for _ in range(0, N)]
    upper = [-A + B for _ in range(0, N)]

    data = lower, main, upper
    diags = np.array([-1, 0, 1])
    matrix = sparse.spdiags(data, diags, N, N).todense()

    # Set values for periodic boundary conditions
    matrix[0, N - 1] = -A - B
    matrix[N - 1, 0] = -A + B

    # Initialize RHS
    RHS = np.array(Phi_old)

    t = 0
    while t <= tau:
        # Enforce our periodic boundary condition
        RHS[0] = (A * (Phi_old[1] - 2 * Phi_old[0] + Phi_old[-1]) -
                  B * (Phi_old[1] - Phi_old[-1]) +
                  Phi_old[0])

        for i in range(1, N - 1):
            RHS[i] = (A * (Phi_old[i + 1] - 2 * Phi_old[i] + Phi_old[i - 1]) -
                      B * (Phi_old[i + 1] - Phi_old[i - 1]) +
                      Phi_old[i])

        # Enforce our periodic boundary condition
        RHS[-1] = (A * (Phi_old[0] - 2 * Phi_old[-1] + Phi_old[-2]) -
                   B * (Phi_old[0] - Phi_old[-2]) +
                   Phi_old[-1])

        # Solve matrix
        Phi = np.linalg.solve(matrix, RHS)

        Phi_old = np.array(Phi)
        t += dt

    return np.array(Phi_old)


def QUICK(Phi, c):
    D, dt, dx, u, tau = c.D, c.dt, c.dx, c.u, c.tau

    N = len(Phi)
    Phi = np.array(Phi)
    Phi_old = np.array(Phi)

    A = D * dt / dx**2
    B = u * dt / (8 * dx)

    t = 0
    while t <= tau:
        for i in range(2, N - 1):
            Phi[i] = (A * (Phi_old[i + 1] - 2 * Phi_old[i] + Phi_old[i - 1]) -
                      B * (3 * Phi_old[i + 1] + Phi_old[i - 2] - 7 * Phi_old[i - 1] + 3 * Phi_old[i]) +
                      Phi_old[i])

        Phi[-1] = (A * (Phi_old[1] - 2 * Phi_old[-1] + Phi_old[-2]) -
                   B * (3 * Phi_old[1] + Phi_old[-3] - 7 * Phi_old[-2] + 3 * Phi_old[-1]) +
                   Phi_old[-1])
        Phi[0] = Phi[-1]
        Phi[1] = (A * (Phi_old[2] - 2 * Phi_old[1] + Phi_old[0]) -
                  B * (3 * Phi_old[2] + Phi_old[-2] - 7 * Phi_old[0] + 3 * Phi_old[1]) +
                  Phi_old[1])

        # Increment
        Phi_old = np.array(Phi)
        t += dt

    return np.array(Phi_old)


def stability(c):
    C, s, u = c.C, c.s, c.u

    FTCS = C <= np.sqrt(2 * s * u) and s <= 0.5
    Upwind = C + 2*s <= 1
    Trapezoidal = True
    QUICK = C <= min(2 - 4 * s, np.sqrt(2 * s))

    # print('C = ', C, ' s = ', s)
    # print('FTCS: ' + str(FTCS))
    # print('Upwind: ' + str(Upwind))
    # print('Trapezoidal: ' + str(Trapezoidal))
    # print('QUICK: ' + str(QUICK))

    return [FTCS, Upwind, Trapezoidal, QUICK]


def calc_stability(C, s, solver):
    results = []
    for C_i, s_i in zip(C, s):
        out = generate_solutions(C_i, s_i, find_order=True)
        results.append(out)

    # Sort and convert
    results.sort(key=lambda x: x[0])
    results = np.array(results)

    # Pull out data
    x = results[:, 0]
    t = results[:, 1]
    RMS_FTCS = results[:, 2]
    RMS_Upwind = results[:, 3]
    RMS_Trapezoidal = results[:, 4]
    RMS_QUICK = results[:, 5]

    # Plot effective orders
    rms_list = [(RMS_FTCS, 'FTCS'),
                (RMS_Upwind, 'Upwind'),
                (RMS_Trapezoidal, 'Trapezoidal'),
                (RMS_QUICK, 'QUICK')]

    for rms in rms_list:
        if rms[1] == solver:
            plot_order(x, t, rms)


def generate_solutions(C, s, find_order=False):
    c = Config(C, s)

    # Spit out some stability information
    stable = stability(c)

    # Initial Condition with boundary conditions
    Phi_initial = np.sin(c.k * c.x)

    # Analytic Solution
    Phi_analytic = Analytic(c)

    # Explicit Solution
    Phi_ftcs = FTCS(Phi_initial, c)

    # Upwind Solution
    Phi_upwind = Upwind(Phi_initial, c)

    # Trapezoidal Solution
    Phi_trapezoidal = Trapezoidal(Phi_initial, c)

    # QUICK Solution
    Phi_quick = QUICK(Phi_initial, c)

    # Save group comparison
    solutions = [(Phi_ftcs, 'FTCS'),
                 (Phi_upwind, 'Upwind'),
                 (Phi_trapezoidal, 'Trapezoidal'),
                 (Phi_quick, 'QUICK')]

    if not find_order:
        # Save individual comparisons
        save_figure(c.x, Phi_analytic, Phi_ftcs,
                    'FTCS ' + str(C) + ' ' + str(s), stable[0])
        save_figure(c.x, Phi_analytic, Phi_upwind,
                    'Upwind ' + str(C) + ' ' + str(s), stable[1])
        save_figure(c.x, Phi_analytic, Phi_trapezoidal,
                    'Trapezoidal ' + str(C) + ' ' + str(s), stable[2])
        save_figure(c.x, Phi_analytic, Phi_quick,
                    'QUICK ' + str(C) + ' ' + str(s), stable[3])

        # and group comparisons
        save_state(c.x, Phi_analytic, solutions, str(C) + ' ' + str(s))
        save_state_error(c.x, Phi_analytic, solutions, str(C) + ' ' + str(s))

    NRMS = []
    for solution in solutions:
        err = solution[0] - Phi_analytic
        NRMS.append(np.sqrt(np.mean(np.square(err)))/(max(Phi_analytic) - min(Phi_analytic)))

    return [c.dx, c.dt, NRMS[0], NRMS[1], NRMS[2], NRMS[3]]


def main():
    # Cases
    C = [0.1,   0.5,   2, 0.5, 0.5]
    s = [0.25, 0.25, .25, 0.5,   1]
    for C_i, s_i in zip(C, s):
        generate_solutions(C_i, s_i)

    # Stable values for each case to find effective order of methods
    C = [0.10, 0.50, 0.40, 0.35, 0.5]
    s = [0.25, 0.25, 0.25, 0.40, 0.5]
    calc_stability(C, s, 'FTCS')

    C = [0.1, 0.2, 0.3, 0.05, 0.1]
    s = [0.4, 0.3, 0.2,  0.1, 0.1]
    calc_stability(C, s, 'Upwind')

    C = [0.5, 0.6, 0.7, 0.8, 0.9]
    s = [0.25, 0.25, 0.25, 0.25, 0.25]
    calc_stability(C, s, 'Trapezoidal')

    C = [0.25, 0.4, 0.5, 0.6,  0.7]
    s = [0.25, 0.25, 0.25, 0.25, 0.25]
    calc_stability(C, s, 'QUICK')


if __name__ == "__main__":
    main()
