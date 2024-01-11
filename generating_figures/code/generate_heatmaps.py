# AUTHOR: GASTON LUCA CASTILLO MORO
# DATE: JANUARY 11, 2024

import matplotlib.pyplot as plt
import matplotlib.colors as colors

plt.rcParams['text.usetex'] = True # this allows the use of LaTeX for the labels. Further commands may be needed depending on Python environment and you should refer to the appropriate documentation if this is the case

def heatmap_plotting_1v1_comparison(data, p_range, x_ticks, y_ticks, end_pd = 120, end_pe = 120, colormap = 'magma', dividend_name = 'WSP', divisor_name = 'YFI'):
  # the data argument was the final probabilities that were calculated with the data-generating code
  X, Y = np.meshgrid(p_range[:end_pd], p_range[:end_pd])

  # log2 transform the data for better determination of patterns
  data = np.log2(data)[:end_pd, :end_pe]

  # Create heatmap with explicit coordinates
  plt.pcolormesh(X, Y, data, vmin = -1.2, vmax = 1.2, cmap=colormap)
  cbar = plt.colorbar()

  # Set log scale for both x and y axes
  plt.xscale('log')
  plt.yscale('log')

  # Add labels
  plt.xlabel('Probability of enabling mutation', labelpad=8, size = 12)
  plt.ylabel('Probability of disabling mutation', labelpad=8, size = 12)

  # change the colorbar
  cbar.set_label(r'$\log_2$ ratio of probability of ' + str(dividend_name) + '/' + str(divisor_name), rotation=90, labelpad=8, size = 12)

  name = str(dividend_name) + ' divided by ' + str(divisor_name) + '.jpg'
  plt.savefig(name)
  plt.show()

def generate_all_1v1_heatmaps():
  yfi = castillom_yfi()
  yfi_probabilities, p_range = retrieve_final_probabilities(yfi.name, 120)
  wsp = castillom_wsp()
  wsp_probabilities, p_range2 = retrieve_final_probabilities(wsp.name, 120)
  sia = castillom_sia()
  sia_probabilities, p_range3 = retrieve_final_probabilities(sia.name, 120)
  morA = castillom_morA()
  morA_probabilities, p_range4 = retrieve_final_probabilities(morA.name, 120)

  # for simplicity, we assumed that pd_range and pe_range were the same and that all the p_range's were the same (as that was what we used when generating the figures)

  heatmap_plotting_1v1_comparison(np.divide(wsp_probabilities, yfi_probabilities), p_range, [10e-7, 10e-6, 10e-5, 10e-4, 10e-3, 10e-2], [10e-7, 10e-6, 10e-5, 10e-4, 10e-3, 10e-2], 120, 120, colormap = 'RdBu', dividend_name = 'Wsp', divisor_name = 'Yfi')
  heatmap_plotting_1v1_comparison(np.divide(wsp_probabilities, sia_probabilities), p_range, [10e-7, 10e-6, 10e-5, 10e-4, 10e-3, 10e-2], [10e-7, 10e-6, 10e-5, 10e-4, 10e-3, 10e-2], 120, 120, colormap = 'RdBu', dividend_name = 'Wsp', divisor_name = 'Sia')
  heatmap_plotting_1v1_comparison(np.divide(wsp_probabilities, morA_probabilities), p_range, [10e-7, 10e-6, 10e-5, 10e-4, 10e-3, 10e-2], [10e-7, 10e-6, 10e-5, 10e-4, 10e-3, 10e-2], 120, 120, colormap = 'RdBu', dividend_name = 'Wsp', divisor_name = 'MorA')
  heatmap_plotting_1v1_comparison(np.divide(yfi_probabilities, sia_probabilities), p_range, [10e-7, 10e-6, 10e-5, 10e-4, 10e-3, 10e-2], [10e-7, 10e-6, 10e-5, 10e-4, 10e-3, 10e-2], 120, 120, colormap = 'RdBu', dividend_name = 'Yfi', divisor_name = 'Sia')
  heatmap_plotting_1v1_comparison(np.divide(yfi_probabilities, morA_probabilities), p_range, [10e-7, 10e-6, 10e-5, 10e-4, 10e-3, 10e-2], [10e-7, 10e-6, 10e-5, 10e-4, 10e-3, 10e-2], 120, 120, colormap = 'RdBu', dividend_name = 'Yfi', divisor_name = 'MorA')
  heatmap_plotting_1v1_comparison(np.divide(sia_probabilities, morA_probabilities), p_range, [10e-7, 10e-6, 10e-5, 10e-4, 10e-3, 10e-2], [10e-7, 10e-6, 10e-5, 10e-4, 10e-3, 10e-2], 120, 120, colormap = 'RdBu', dividend_name = 'Sia', divisor_name = 'MorA')

if name == '__main__':
  generate_all_1v1_heatmaps()
