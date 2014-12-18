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
fig_height_in = fig_width_in * golden_ratio   # figure height in inches
fig_dims = [fig_width_in, fig_height_in]  # fig dims as a list

arr_befor = []
arr_after = []
for i in range(201):
    j = str(i)

    # Load data
    path = 'M=100/Run0.5-100.0/postProcessing/sets/' + j + '/'
    before = np.loadtxt(path + 'before_Ux_Uy.xy', delimiter=' ', unpack=True)
    after  = np.loadtxt(path + 'after_Ux_Uy.xy',  delimiter=' ', unpack=True)

    Y, UX, UY = 0, 1, 2
    before_Y,  after_Y  = before[Y],  after[Y]
    before_UX, after_UX = before[UX], after[UX]
    # before_UY, after_UY = before[UY], after[UY]

    # Final flow should be linear
    linear = np.linspace(0, 1, len(before_UX))

    # Calculate RMS
    err_befor, err_after = before_UX - linear, after_UX - linear
    arr_befor.append(np.sqrt(np.mean(np.square(err_befor))))
    arr_after.append(np.sqrt(np.mean(np.square(err_after))))


# Plot
plt.figure(figsize=fig_dims)
plt.plot(arr_befor, 'b+', label='before')
plt.plot(arr_after, 'rx', label='after')
plt.yscale('log')
plt.legend()
plt.xlabel('Time (s)')
plt.ylabel('NRMS')

# plt.figure(figsize=fig_dims)
# plt.subplot(2, 1, 1)

# plt.plot(befor_error, before_Y, 'b+', label='before')
# plt.plot(after_error,  after_Y, 'rx', label='after')
# # plt.plot(middle_Y, middle_UX,  'g:', label='middle')
# plt.ylabel('y (m)')
# plt.xlabel('Error in Ux (m/s)')
# # plt.xlim([0, 1])
# plt.legend(loc='upper left')

# plt.subplot(2, 1, 2)
# plt.plot(before_Y, before_UY, 'b+',  label='before')
# plt.plot(after_Y,  after_UY,  'rx', label='after')
# # plt.plot(middle_Y, middle_UY,  'g:', label='middle')
# plt.xlabel('Y')
# plt.ylabel('Uy')
# plt.xlim([0, 1])
# plt.legend(loc='upper left')
plt.show()
