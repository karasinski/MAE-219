import numpy as np
import matplotlib.pyplot as plt

# Problem Parameters
L = 1.            # Domain lenghth       [n.d.]
T0 = 0.           # Initial temperature  [n.d.]
T1 = 1.           # Boundary temperature [n.d.]
t_start = 0.
t_end = 0.1
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

plt.axis([0, L, T0, T1])
plt.xlabel('Length [nd]')
plt.ylabel('Temperature [nd]')

while time <= t_end:
    for i in range(1, N - 1):
        Tnew[i] = s * Told[i + 1] + (1 - 2.0 * s) * Told[i] + s * Told[i - 1]

    plt.plot(x, Tnew, linewidth=1)
    time = time + dt
    Told = Tnew

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

a = [0.0,  2.0, 1.0]  # first number must be zero
b = [1.0,  1.0, 1.0]
c = [2.0,  1.0, 0.0]  # last number must be zero
d = [0.0,  1.0, 0.0]

A = np.array([[b[0], c[0],   0],
              [a[1], b[1], c[1]],
              [  0,  a[2], b[2]]])

# Python's linalg solver
x = np.linalg.solve(A, d)

# My solver
solved = TDMAsolver(a, b, c, d)

# Does my solution solve the problem?
print(np.allclose(np.dot(A, solved), d))

# Does my solution match the linalg solver?
print(np.allclose(x, solved))
