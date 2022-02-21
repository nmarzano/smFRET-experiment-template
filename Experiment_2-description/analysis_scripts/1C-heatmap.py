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
import glob as glob

from Utilities.Data_analysis import file_reader

input_folder = 'Figure_1-Validation-of-Fluc-sensor/raw_data/200715_Fluc_FoldandUnfoldtraces'
output_folder = 'Figure_1-Validation-of-Fluc-sensor/python_results'
exposure = 0.500  # in seconds
frame_rate = 1/exposure

data_paths = {
    'Heatmap':('Heatmap', 'Figure_1-Validation-of-Fluc-sensor/raw_data/200715_Fluc_FoldandUnfoldtraces')
}
initial = file_reader(input_folder, 'heatmap', frame_rate)


def remove_outliers(compiled, plot_type, data_type = "raw"):
    """[removes outliers from dataframe]

    Args:
        compiled ([dataframe]): [raw dataframe containing outliers to be removed]
        plot_type ([str]): [string can either be 'hist' for histogram data or 'TDP' for TDP data]
        data_type (str, optional): [removes either raw FRET values or 'idealized' FRET values]. Defaults to "raw".

    Returns:
        [dataframe]: [returns cleaned data without outliers]
    """
    if plot_type == 'hist':
        if data_type == "raw":
            rawFRET = compiled[(compiled['FRET'] > -0.5) & (compiled['FRET'] < 1.5)].copy()
            return rawFRET
        if data_type == "idealized":
            idealizedFRET = compiled[(compiled[4] > -0.5) & (compiled[4] < 1.5)].copy()
            return idealizedFRET
    elif plot_type == 'TDP':
        outliers = compiled[(compiled["FRET before transition"] < -0.5)|(compiled["FRET before transition"] > 1.5)|(compiled["FRET after transition"] < -0.5) | (compiled["FRET after transition"] > 1.5)].index
        compiled.drop(outliers, inplace = True)
        return compiled
    else:
        print('invalid plot type, please set plot_type as "hist" or "TDP" - you idiot')

compiled_df = []
for data_name, (label, data_path) in data_paths.items():
    imported_data = file_reader(input_folder, 'heatmap', frame_rate)
    cleaned_raw = remove_outliers(imported_data, 'hist')    #### add "idealized" after imported_data to get idealized histograms
    cleaned_raw["treatment_name"] = data_name
    compiled_df.append(cleaned_raw)
compiled_df = pd.concat(compiled_df)   #### .rename(columns = {1:"test", 3:"test2"}) ## can rename individually if needed
compiled_df.columns = ["frames", "donor", "acceptor", "FRET", "idealized FRET", 'molecule_number', 'time', "treatment_name"]




plt.rcParams['svg.fonttype'] = 'none'
sns.set(style="whitegrid")
g = sns.JointGrid(data = initial, x='time', y='FRET', xlim = (0,200), ylim = (0, 1.1))
g.plot_joint(plt.hexbin, gridsize=(50,50), cmap='ocean_r', mincnt=0, bins=80)
g.plot_marginals(sns.histplot, kde=True, bins=20)
identifier = input_folder.split('/')[-1]
plt.savefig(f'{output_folder}/Heatmap_{identifier}.svg', dpi = 600)


# # creates the grey version of the colormap,  each R, G and B can be adjusted to change the colour
# N = 256
# yellow = np.ones((N, 4))
# yellow[:, 0] = np.linspace(64/256, 1, N) # R = 255
# yellow[:, 1] = np.linspace(64/256, 1, N) # G = 232
# yellow[:, 2] = np.linspace(64/256, 1, N)  # B = 11
# yellow_cmp = ListedColormap(yellow)

# # creates the blue version of the colormap,  each R, G and B can be adjusted to change the colour
# red = np.ones((N, 4))
# red[:, 0] = np.linspace(0/256, 1, N)    # R = 255
# red[:, 1] = np.linspace(139/256, 1, N)  # G = 232
# red[:, 2] = np.linspace(226/256, 1, N)  # B = 11
# red_cmp = ListedColormap(red)
# # stacks the two colormaps that have just been defined above. Changing the last number in the grey and blue sequences will change 
# # the relative proportion of the concatenated colormap. The sum must equal 258.
# newcolors2 = np.vstack((yellow_cmp(np.linspace(0, 1, 32)),
#                        red_cmp(np.linspace(1, 0, 226))))

# # Renames the new concatenated colormap to a cmap
# double = ListedColormap(newcolors2, name='double')

# # convert to matplotlib colormap
# # cmap = mpl.colors.ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])

# gs = gridspec.GridSpec(1, 2)
# cbaxes = inset_axes(g, width="30%", height="3%", loc=3) 
# plt.colorbar(cax=cbaxes, ticks=[0.,1], orientation='horizontal')
# plt.show()
# g.ax_joint.set_facecolor("#404040") # to change the colour of the background
