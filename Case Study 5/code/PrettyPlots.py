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
    plt.ylim([1E4, 1E8])
    plt.yticks([1E4, 2E4, 3E4, 4E4, 5E4,
                1E5, 2E5, 3E5, 4E5, 5E5,
                1E6, 2E6, 3E6, 4E6, 5E6,
                1E7, 2E7, 3E7, 4E7, 5E7, 1E8],
               ['E+4', '2', '3', '4', '5',
                'E+5', '2', '3', '4', '5',
                'E+6', '2', '3', '4', '5',
                'E+7', '2', '3', '4', '5', 'E+8'])
    plt.legend(loc='lower right')

    plt.text(30.5, 1.5e+4, 't=2', fontsize=9, family='serif')
    plt.text(30.5, 3.5e+5, 't=0', fontsize=9, family='serif')
    plt.text(30.5, 1.5e+6, 't=9', fontsize=9, family='serif')
    plt.text(30.5, 1.e+7,  't=4', fontsize=9, family='serif')
    plt.text(30.5, 2.8e+7, 't=7', fontsize=9, family='serif')
    plt.text(30.5, 6.e+7,  't=6', fontsize=9, family='serif')

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
    plt.yscale('log')
    plt.ylim([4.E11, 1.2E12])
    plt.yticks([4E11, 5E11, 6E11, 7E11, 8E11, 9E11, 1E12, 1.2E12],
               ['4.00e+11', '5', '6', '7', '8', '9', 'E+12', '1.20e+12'])
    plt.legend(loc='best')

    # Left side text
    plt.text(28.25, 4.65e+11, 't=0,2,4', fontsize=9, family='serif')
    plt.text(29, 5.4e+11,     't=6',     fontsize=9, family='serif')
    plt.text(29, 5.7e+11,     't=7',     fontsize=9, family='serif')

    # Right side text
    plt.text(50.25, 4.95e+11, 't=0',  fontsize=9, family='serif')
    plt.text(50.25, 5.35e+11, 't=2',  fontsize=9, family='serif')
    plt.text(50.25, 5.6e+11,  't=4',  fontsize=9, family='serif')
    plt.text(50.25, 6.1e+11,  't=6',  fontsize=9, family='serif')
    plt.text(50.25, 6.4e+11,  't=7',  fontsize=9, family='serif')
    plt.text(50.25, 6.7e+11,  't=9',  fontsize=9, family='serif')
    plt.text(50.25, 7.0e+11,  't=12', fontsize=9, family='serif')
    plt.text(50.25, 7.3e+11,  't=18', fontsize=9, family='serif')
    plt.text(50.25, 7.6e+11,  't=24', fontsize=9, family='serif')

    save_name = integrator + ' c2.pdf'
    save_plot(save_name)


def plot_40km(t, c1_40km, c2_40km, integrator):
    c2_40km_scaled = [1E-4 * val for val in c2_40km]
    days = [val / 86400. for val in t]

    plt.figure(figsize=fig_dims)
    plt.plot(days, c1_40km, label='$c_1$')
    plt.plot(days, c2_40km_scaled, label='$c_2$ * 1E-4')
    plt.xlabel('t (days)')
    plt.yscale('log')
    plt.ylim([2.E6, 2E8])
    plt.yticks([2E6, 3E6, 4E6, 5E6, 1E7, 2E7, 3E7, 4E7, 5E7, 1E8, 2E8],
               ['2', '3', '4', '5', 'E+7', '2', '3', '4', '5', 'E+8', '2'])
    plt.xlim([0, days[-1]])
    plt.legend(loc='lower right')
    save_name = integrator + ' time.pdf'
    save_plot(save_name)
