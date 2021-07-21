from scipy.stats import norm
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import pandas as pd
import numpy as np
import Utilities.Data_analysis as uda

output_folder = "Experiment_X-description/python_results"
FRET_thresh = 0.5 #### FRET value at which to filter data above or below. IF CHANGED, WILL NEED TO CHANGE ALL 0.5 VALUES (E.G. BELOW IN HEADERS) TO THE NEW VALUE
############# NEED TO MODIFY BASED ON DATASET 
data_paths_violin = {
    "treatment_label_1":"data_directory_1/Filtered_dwelltime_treatment_1.csv",
    "treatment_label_2":"data_directory_2/Filtered_dwelltime_treatment_2.csv",
    "treatment_label_3":"data_directory_3/Filtered_dwelltime_treatment_3.csv",
    "treatment_label_4":"data_directory_4/Filtered_dwelltime_treatment_4.csv",
    "treatment_label_5":'data_directory_5/Filtered_dwelltime_treatment_5.csv'
}


def compile(df, data_name):
    violin_data_lowtolow = pd.DataFrame(df[f"< {FRET_thresh} to < {FRET_thresh}"])
    violin_data_lowtolow.columns = ["y_axis"]
    violin_data_lowtolow["transition_type"] = f"< {FRET_thresh} to < {FRET_thresh}"
    violin_data_lowtolow["treatment"] = data_name

    violin_data_lowtohigh = pd.DataFrame(df[f"< {FRET_thresh} to > {FRET_thresh}"])
    violin_data_lowtohigh.columns = ["y_axis"]
    violin_data_lowtohigh["transition_type"] = f"< {FRET_thresh} to > {FRET_thresh}"
    violin_data_lowtohigh["treatment"] = data_name

    violin_data_hightohigh = pd.DataFrame(df[f"> {FRET_thresh} to > {FRET_thresh}"])
    violin_data_hightohigh.columns = ["y_axis"]
    violin_data_hightohigh["transition_type"] = f"> {FRET_thresh} to > {FRET_thresh}"
    violin_data_hightohigh["treatment"] = data_name

    violin_data_hightolow = pd.DataFrame(df[f"> {FRET_thresh} to < {FRET_thresh}"])
    violin_data_hightolow.columns = ["y_axis"]
    violin_data_hightolow["transition_type"] = f"> {FRET_thresh} to < {FRET_thresh}"
    violin_data_hightolow["treatment"] = data_name
    merged = pd.concat([violin_data_lowtolow, violin_data_lowtohigh, violin_data_hightohigh, violin_data_hightolow])
    return merged


test = []
for data_name, data_path in data_paths_violin.items():
    data = uda.file_reader(data_path, 'other')
    compiled_data = compile(data, data_name)
    test.append(compiled_data)
final = pd.concat(test)
final["y_axis_log10"] = np.log10(final['y_axis']) ## if need to plot in log scale

############# Code to plot data
############# NEED TO MODIFY BASED ON DATASET but can COPY FROM 1A-PLOT-HISTOGRAM
colors_violin = {
    "Native": "black", 
    "Spontaneous":"darkorange", 
    "treatment_3":"skyblue", 
    "treatment_4":"royalblue", 
    "treatment_5":"darkviolet"
}


def plot_violin(data, scale = "y_axis"):
    if scale == "y_axis":
        plt.rcParams['svg.fonttype'] = 'none'
        plot1 = plt.figure(figsize = (12, 6))
        sns.set(style = "darkgrid", font_scale = 1.5)
        sns.violinplot(
            data=data, 
            x= "transition_type", 
            y = "y_axis",
            palette= colors_violin, 
            hue = "treatment", 
            log_scale = True,
            cut = 0)
        plt.ylabel("Residence time (s)")
        plt.xlabel("Transition class")
        plt.legend(loc = "upper right")
        plot1.savefig(f"{output_folder}/Violin_plot_normal.svg", dpi = 600)
        plt.show()
    if scale == "y_axis_log10":
        plt.rcParams['svg.fonttype'] = 'none'
        plot2 = plt.figure(figsize = (12, 6))
        sns.set(style = "darkgrid", font_scale = 1.5)
        sns.violinplot(
            data=data, 
            x= "transition_type", 
            y = "y_axis_log10",
            palette= colors_violin, 
            hue = "treatment", 
            log_scale = True)
        plt.ylabel("Log residence time (s)")
        plt.xlabel("Transition class")
        plt.legend(loc = "upper left", bbox_to_anchor = (1,1), ncol =1)
        plot2.savefig(f"{output_folder}/Violin_plot_log.svg", dpi = 600)
        plt.show()
    if scale == 'split':     
        f, (ax_top, ax_bottom) = plt.subplots(ncols=1, nrows=2, sharex=True, gridspec_kw={'hspace':0.05})
        sns.violinplot(x="transition_type", y="y_axis", hue="treatment",data=final, ax=ax_top, palette = colors_violin)
        sns.violinplot(x="transition_type", y="y_axis", hue="treatment",data=final, ax=ax_bottom, cut = 0, palette = colors_violin)
        ax_top.set_ylim(bottom=40)   # those limits are fake
        ax_bottom.set_ylim(0,40)
        sns.despine(ax=ax_bottom)
        sns.despine(ax=ax_top, bottom=True)
        ax = ax_top
        d = .015  # how big to make the diagonal lines in axes coordinates
        # arguments to pass to plot, just so we don't keep repeating them
        kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
        ax.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
        ax2 = ax_bottom
        kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
        ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
        #remove one of the legend
        ax_bottom.legend_.remove()
        f.savefig(f"{output_folder}/Violin_plot_splitaxis.svg", dpi = 600)
        plt.show()

plot_violin(final, 'split')


############# Code to generate and collate summary statistics of dwell times
mean = final.groupby(['treatment', 'transition_type']).mean()
sem =  final.groupby(['treatment', 'transition_type']).sem()
N = final.groupby(['treatment', 'transition_type']).count()
collated = pd.concat([mean,sem, N], axis = 1)
collated.drop([col for col in collated.columns.tolist() if 'y_axis_log10' in col], axis = 1, inplace = True)
collated.columns = ['mean_residence_time', 'sem', 'n']
collated.reset_index(inplace = True)
collated.to_csv(f"{output_folder}/summary.csv", index = False)


##### "Pastel1" - this is quite a nice colour scheme