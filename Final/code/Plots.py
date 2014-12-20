import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

# Configure figures for production
WIDTH = 495.0  # the number latex spits out
FACTOR = 1.0   # the fraction of the width the figure should occupy
fig_width_pt = WIDTH * FACTOR

inches_per_pt = 1.0 / 72.27
golden_ratio = (np.sqrt(5) - 1.0) / 2.0      # because it looks good
fig_width_in = fig_width_pt * inches_per_pt  # figure width in inches
fig_height_in = fig_width_in * golden_ratio  # figure height in inches
fig_dims = [fig_width_in, fig_height_in]     # fig dims as a list

ARs = ["0.5", "2.0", "5.0"]
colors = ['r', 'g', 'b']
plt.figure(figsize=fig_dims)
for AR, color in zip(ARs, colors):
    arr_befor = []
    arr_after = []
    for i in range(150):
        # Load data
        path = 'M=100/Run' + AR + '-100.0/postProcessing/sets/' + str(i) + '/'
        before = np.loadtxt(path + 'before_Ux_Uy.xy', delimiter=' ', unpack=True)
        after  = np.loadtxt(path + 'after_Ux_Uy.xy',  delimiter=' ', unpack=True)

        Y, UX = 0, 1
        before_Y,  after_Y  = before[Y],  after[Y]
        before_UX, after_UX = before[UX], after[UX]

        # Final flow should be linear
        linear = np.linspace(0, 1, len(before_UX))

        # Calculate RMS
        err_befor, err_after = before_UX - linear, after_UX - linear
        arr_befor.append(np.sqrt(np.mean(np.square(err_befor))))
        arr_after.append(np.sqrt(np.mean(np.square(err_after))))

    # Plot
    plt.plot(arr_befor, color + '+', label='AR=' + AR + ', before')
    plt.plot(arr_after, color + 'x', label='AR=' + AR + ', after')
    plt.yscale('log')
    plt.legend()
    plt.xlim([0, 150])
    plt.xlabel('Time (s)')
    plt.ylabel('NRMS')

plt.show()
plt.close()

path = 'M=100/Run0.5-100.0/postProcessing/sets/30/'
before = np.loadtxt(path + 'before_Ux_Uy.xy', delimiter=' ', unpack=True)
after  = np.loadtxt(path + 'after_Ux_Uy.xy',  delimiter=' ', unpack=True)

Y, UX = 0, 1
before_Y,  after_Y  = before[Y],  after[Y]
before_UX, after_UX = before[UX], after[UX]

plt.figure(figsize=fig_dims)
plt.plot(before_UX, before_Y, 'b+', label='before')
plt.plot( after_UX,  after_Y, 'rx', label='after')
plt.ylabel('y (m)')
plt.xlabel('Error in Ux (m/s)')
plt.legend(loc='upper left')
plt.show()
