def Homework1():
    import numpy as np
    import matplotlib.pyplot as plt

    t = 0  # s
    dt = .1  # s
    t_f = 10  # s

    length = 1  # m
    number_of_points = 10
    dx = length / number_of_points  # m
    x = np.linspace(0, length, number_of_points)

    alpha = 2E-2  # m^2/s
    C = (alpha * dt / dx ** 2)  # better to calculate once than every timestep

    # set up our Temperature list with initial conditions
    Tn_ = np.zeros(number_of_points)
    Tn_[0] = 1
    Tn_[-1] = 1

    history = Tn_
    history = np.vstack((history, Tn_))
    # print(t)

    while t < t_f - dt:
        # make a copy
        Tn = Tn_

        # update all the interior nodes
        for i in range(1, number_of_points - 1):
            Tn[i] = Tn_[i] + C * (Tn_[i + 1] - 2 * Tn_[i] + Tn_[i - 1])

        # update our array and our current time
        Tn_ = Tn
        t += dt

        # append to our history every second (made hellish by floating point)
        if (abs(int(t) - t) < dt):
            history = np.vstack((history, Tn))
            # print(t)
    history = np.vstack((history, Tn))
    # print(t)

    # plot the rod every second
    plt.figure()
    plt.hold(True)
    for i, T in enumerate(history):
        if i == 0:
            continue
        plt.plot(x, history[i], '-', label=str(i - 1) + ' seconds')
    plt.legend(loc='best')
    plt.show()

Homework1()
