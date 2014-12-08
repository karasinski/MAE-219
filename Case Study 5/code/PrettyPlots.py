import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import os

# Configure figures for production
WIDTH = 495.0  # the number latex spits out
FACTOR = 1.0   # the fraction of the width the figure should occupy
fig_width_pt = WIDTH * FACTOR

inches_per_pt = 1.0 / 72.27
golden_ratio = (np.sqrt(5) - 1.0) / 2.0      # because it looks good
fig_width_in = fig_width_pt * inches_per_pt  # figure width in inches
fig_height_in = fig_width_in * golden_ratio   # figure height in inches
fig_dims = [fig_width_in, fig_height_in]  # fig dims as a list


def save_plot(save_name):
    # Save plots
    try:
        os.mkdir('figures')
    except Exception:
        pass

    plt.savefig('figures/' + save_name, bbox_inches='tight')
    plt.close()


def plot_c1(z, initial, c1, labels, integrator):
    plt.figure(figsize=fig_dims)
    plt.plot(z, initial[0::2], label='0 hours')
    for solution, label in zip(c1, labels):
        if "12" in label:
            break
        plt.plot(z, solution, label=label)
    plt.ylabel('$c_1$')
    plt.xlabel('z (km)')
    plt.yscale('log')
    plt.legend()
    save_name = integrator + ' c1.pdf'
    save_plot(save_name)


def plot_c2(z, initial, c2, labels, integrator):
    plt.figure(figsize=fig_dims)
    plt.plot(z, initial[1::2], label='0 hours')
    for solution, label in zip(c2, labels):
        if "240" in label:
            break
        plt.plot(z, solution, label=label)
    plt.ylabel('$c_2$')
    plt.xlabel('z (km)')
    plt.legend()
    save_name = integrator + ' c2.pdf'
    save_plot(save_name)


def plot_40km(t, c1_40km, c2_40km, integrator):
    c2_40km_scaled = [1E-4 * val for val in c2_40km]

    plt.figure(figsize=fig_dims)
    plt.plot(t, c1_40km, label='$c_1$')
    plt.plot(t, c2_40km_scaled, label='$c_2$ * 1E-4')
    plt.xlabel('t (s)')
    plt.legend()
    plt.xlim([0, t[-1]])

    save_name = integrator + ' time.pdf'
    save_plot(save_name)
