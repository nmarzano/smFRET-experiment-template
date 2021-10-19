from scipy.stats import norm
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import pandas as pd
import numpy as np
from matplotlib import rc
import matplotlib as mpl
from matplotlib import cm 
from matplotlib.colors import ListedColormap,LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


# creates the grey version of the colormap,  each R, G and B can be adjusted to change the colour
N = 256
yellow = np.ones((N, 4))
yellow[:, 0] = np.linspace(64/256, 1, N) # R = 255
yellow[:, 1] = np.linspace(64/256, 1, N) # G = 232
yellow[:, 2] = np.linspace(64/256, 1, N)  # B = 11
yellow_cmp = ListedColormap(yellow)

# creates the blue version of the colormap,  each R, G and B can be adjusted to change the colour
red = np.ones((N, 4))
red[:, 0] = np.linspace(0/256, 1, N)    # R = 255
red[:, 1] = np.linspace(139/256, 1, N)  # G = 232
red[:, 2] = np.linspace(226/256, 1, N)  # B = 11
red_cmp = ListedColormap(red)
# stacks the two colormaps that have just been defined above. Changing the last number in the grey and blue sequences will change 
# the relative proportion of the concatenated colormap. The sum must equal 258.
newcolors2 = np.vstack((yellow_cmp(np.linspace(0, 1, 32)),
                       red_cmp(np.linspace(1, 0, 226))))

# Renames the new concatenated colormap to a cmap
double = ListedColormap(newcolors2, name='double')

# convert to matplotlib colormap
# cmap = mpl.colors.ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])

filename = 'Experiment_1-High-conc-chap_removal-with-GrpE/raw_data/GrpE_only_tracesHeatmap.dat'

#Read dat files
x = pd.read_fwf(filename)


x.columns = ['Time (s)', 'FRET']


T = x['Time (s)']/5
FRET = x['FRET']


plt.rcParams['svg.fonttype'] = 'none'

g = sns.JointGrid(x=T, y=FRET, xlim = (0,200), ylim = (0, 1.1))
g.plot_joint(plt.hexbin, gridsize=(75,75), cmap=double, mincnt=0, bins=100)
g.plot_marginals(sns.histplot, kde=True, bins=20)
g.ax_joint.set_facecolor("#404040") # to change the colour of the background




# gs = gridspec.GridSpec(1, 2)
# cbaxes = inset_axes(g, width="30%", height="3%", loc=3) 
# plt.colorbar(cax=cbaxes, ticks=[0.,1], orientation='horizontal')
# plt.show()

plt.savefig('Experiment_1-High-conc-chap_removal-with-GrpE/python_results/Heatmap.svg', dpi = 600)