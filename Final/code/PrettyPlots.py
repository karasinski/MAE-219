# import numpy as np
# import matplotlib.pyplot as plt
# import os


# # Configure figures for production
# WIDTH = 495.0  # width of one column
# FACTOR = 1.0   # the fraction of the width the figure should occupy
# fig_width_pt  = WIDTH * FACTOR

# inches_per_pt = 1.0 / 72.27
# golden_ratio  = (np.sqrt(5) - 1.0) / 2.0       # because it looks good
# fig_width_in  = fig_width_pt * inches_per_pt   # figure width in inches
# fig_height_in = fig_width_in * golden_ratio    # figure height in inches
# fig_dims      = [fig_width_in, fig_height_in]  # fig dims as a list


# def sigma_xx(x):
#     return 1E4*(1+(0.125/(x**2))+(0.09375/(x**4)))


# def sigma_yy(x):
#     return 1E4*((0.125/(x**2))-(0.09375/(x**4)))


# def plot_xx(widths, meshes):
#     # Format plot
#     plt.figure(figsize=fig_dims)
#     plt.xlabel('Distance, y (m)')
#     plt.ylabel('Stress ($\sigma_{xx}$)$_{x=0}$(kPa)')
#     title = 'Normal stress along the vertical symmetry'
#     x = np.linspace(0.5, 2)
#     sigmaxx = sigma_xx(x)

#     plt.plot(x, sigmaxx, '-k', label='Analytic Solution')
#     plt.xlim(0.5, 2)

#     for width, mesh in zip(widths, meshes):
#         path = "Run" + str(width) + '-' + str(mesh) + '/postProcessing/sets/100/leftPatch_sigmaxx_sigmayy.xy'
#         data = np.loadtxt(path)

#         if widths.count(widths[0]) == len(widths):
#             label = 'Explicit Solution ($n=' + str(int(mesh)) + '$)'
#         else:
#             label = 'Explicit Solution ($x=' + str(int(2*width)) + '$)'
#         plt.plot(data[:, 0], data[:, 1], '--', markersize=5, label=label)

#     if widths.count(widths[0]) == len(widths):
#         title += ' ($x=' + str(int(2*width)) + '$)'
#     else:
#         title += ' ($n=' + str(int(mesh)) + '$)'

#     plt.title(title)
#     plt.legend(loc='best')

#     # Save plots
#     save_name = 'result-x-' + str(widths) + str(meshes) + '.pdf'
#     try:
#         os.mkdir('figures')
#     except Exception:
#         pass

#     plt.savefig('figures/' + save_name, bbox_inches='tight')
#     plt.clf()


# def plot_xx_err(widths, meshes):
#     # Format plot
#     plt.figure(figsize=fig_dims)
#     plt.xlabel('Distance, y (m)')
#     plt.ylabel('Error in Stress ($\sigma_{xx}$)$_{x=0}$(kPa)')
#     plt.title('Error in Normal stress along the vertical symmetry')
#     plt.xlim(0.5, 2)

#     for width, mesh in zip(widths, meshes):
#         path = "Run" + str(width) + '-' + str(mesh) + '/postProcessing/sets/100/leftPatch_sigmaxx_sigmayy.xy'
#         data = np.loadtxt(path)

#         if widths.count(widths[0]) == len(widths):
#             label = 'Explicit Solution ($n=' + str(int(mesh)) + '$)'
#         else:
#             label = 'Explicit Solution ($x=' + str(int(2*width)) + '$)'

#         x = data[:, 0]
#         sigmaxx = sigma_xx(x)
#         err = data[:, 1] - sigmaxx

#         RMS = np.sqrt(np.mean(np.square(err)))/(max(sigmaxx) - min(sigmaxx))
#         print('x err', width, mesh, '{0:.3e}'.format(RMS))

#         plt.plot(x, err, '--', markersize=5, label=label)

#     plt.legend(loc='best')

#     # Save plots
#     save_name = 'error-x-' + str(widths) + str(meshes) + '.pdf'
#     try:
#         os.mkdir('figures')
#     except Exception:
#         pass

#     plt.savefig('figures/' + save_name, bbox_inches='tight')
#     plt.clf()


# def plot_yy(widths, meshes):
#     # Format plot
#     plt.figure(figsize=fig_dims)
#     plt.xlabel('Distance, x (m)')
#     plt.ylabel('Stress ($\sigma_{yy}$)$_{y=0}$(kPa)')
#     title = 'Normal stress along the horizontal symmetry'
#     y = np.linspace(0.5, 2)
#     sigmayy = sigma_yy(y)

#     plt.plot(y, sigmayy, '-k', label='Analytic Solution')
#     plt.xlim(0.5, 2)

#     for width, mesh in zip(widths, meshes):
#         path = "Run" + str(width) + '-' + str(mesh) + '/postProcessing/sets/100/downPatch_sigmaxx_sigmayy.xy'
#         data = np.loadtxt(path)

#         if widths.count(widths[0]) == len(widths):
#             label = 'Explicit Solution ($n=' + str(int(mesh)) + '$)'
#         else:
#             label = 'Explicit Solution ($x=' + str(int(2*width)) + '$)'
#         plt.plot(data[:, 0], data[:, 2], '--', markersize=5, label=label)

#     if widths.count(widths[0]) == len(widths):
#         title += ' ($x=' + str(int(2*width)) + '$)'
#     else:
#         title += ' ($n=' + str(int(mesh)) + '$)'

#     plt.title(title)
#     plt.legend(loc='best')

#     # Save plots
#     save_name = 'result-y-' + str(widths) + str(meshes) + '.pdf'
#     try:
#         os.mkdir('figures')
#     except Exception:
#         pass

#     plt.savefig('figures/' + save_name, bbox_inches='tight')
#     plt.clf()


# def plot_yy_err(widths, meshes):
#     # Format plot
#     plt.figure(figsize=fig_dims)
#     plt.xlabel('Distance, x (m)')
#     plt.ylabel('Error in Stress ($\sigma_{yy}$)$_{y=0}$(kPa)')
#     plt.title('Error in Normal stress along the horizontal symmetry')
#     plt.xlim(0.5, 2)

#     for width, mesh in zip(widths, meshes):
#         path = "Run" + str(width) + '-' + str(mesh) + '/postProcessing/sets/100/downPatch_sigmaxx_sigmayy.xy'
#         data = np.loadtxt(path)

#         if widths.count(widths[0]) == len(widths):
#             label = 'Explicit Solution ($n=' + str(int(mesh)) + '$)'
#         else:
#             label = 'Explicit Solution ($x=' + str(int(2*width)) + '$)'

#         y = data[:, 0]
#         sigmayy = sigma_yy(y)
#         err = data[:, 2] - sigmayy

#         RMS = np.sqrt(np.mean(np.square(err)))/(max(sigmayy) - min(sigmayy))
#         print('y err', width, mesh, '{0:.3e}'.format(RMS))

#         plt.plot(y, err, '--', markersize=5, label=label)

#     plt.legend(loc='best')

#     # Save plots
#     save_name = 'error-y-' + str(widths) + str(meshes) + '.pdf'
#     try:
#         os.mkdir('figures')
#     except Exception:
#         pass

#     plt.savefig('figures/' + save_name, bbox_inches='tight')
#     plt.clf()


# def generate_plots(widths, meshes):
#     plot_xx(widths, meshes)
#     plot_xx_err(widths, meshes)
#     plot_yy(widths, meshes)
#     plot_yy_err(widths, meshes)

#     print('Plots generated.')
