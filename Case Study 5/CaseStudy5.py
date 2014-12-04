import numpy as np


def gamma(z):
    return 1. - ((z-40.)/10.) ** 2 + 1./2 * ((z-40.)/10.) ** 4


def K(x):
    return 30 + x * dz


def R(n, c, t):
    '''
    Find the reaction rate of the system at state c and time t. It returns R_1
    or R_2, depending on the first input parameter.
    '''

    k_1 = 1.63E-16
    k_2 = 4.66E-16

    a_3 = 22.62
    a_4 = 7.601
    w = np.pi / 43200

    if np.sin(w * t) > 0:
        k_3 = np.exp(-a_3 / np.sin(w * t))
        k_4 = np.exp(-a_4 / np.sin(w * t))
    else:
        k_3 = 0
        k_4 = 0

    if n == 1:
        return (-k_1 * c[0] * c[2]
                - k_2 * c[0] * c[1]
                + 2 * k_3 * c[2]
                + k_4 * c[1])

    elif n == 2:
        return (k_1 * c[0] * c[2]
                - k_2 * c[0] * c[1]
                - k_4 * c[1])
    else:
        raise IndexError('First argument can only be 1 or 2.')


# not sure if these should be c_dot and c or y_dot and y
def solve():
    # print(0, 1)
    c_dot[0] = (dz ** -2 * (K(3/2) * c[2]
                            - (K(3/2) + K(1/2)) * c[0]
                            + K(1/2) * c[2])
                + R(1, c, t))

    c_dot[1] = (dz ** -2 * (K(3/2) * c[3]
                            - (K(3/2) + K(1/2)) * c[1]
                            + K(1/2) * c[3])
                + R(2, c, t))

    for l in range(1, M - 1):
        print(2 * l, 2 * l + 1)

        # the values passed into c on the RHS of these equations is likely wrong
        c_dot[2 * l] = (dz ** -2 * ((K[l + 1./2.] * c[2 * l + 1]
                                    - (K[l + 1./2.] + K[l - 1./2.]) * c[2 * l - 1]
                                    + K[l - 1./2.] * c[2 * l - 3]))
                        + R(1, c[2 * l - 1], c[2 * l], t))

        c_dot[2 * l + 1] = (dz ** -2 * ((K[l + 1./2.] * c[2 * l + 2]
                                        - (K[l + 1./2.] + K[l - 1./2.]) * c[2 * l]
                                        + K[l - 1./2.] * c[2 * l - 2]))
                            + R(2, c[2 * l - 1], c[2 * l], t))

    # print(2*M-2, 2*M-1)


# Basic parameters
M = 50
dz = 20. / M
z = [30. + (j + 1.) * dz for j in range(M)]

# y_1 = 10E6
# y_2 = 10E12
# y_3 = 3.7E16
# y = [y_1, y_2, y_3]

# This generates the initial conditions
c = [0] * (2 * M)
for i in range(M):
    c[2*i] = 10E6 * gamma(30 + (i + 1) * dz)
    c[2*i + 1] = 10E12 * gamma(30 + (i + 1) * dz)

    print(c[2*i], c[2*i + 1])


c_dot = np.zeros(2 * M)
t = 0
