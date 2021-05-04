# library
import seaborn as sns
import pandas as pd
import numpy as np
import glob

FRET_thresh = 0.5

input_folder = 'Experiment_X-description/python_results'
output_folder = "Experiment_X-description/python_results/Heatmaps"

if not os.path.exists(output_folder):
    os.makedirs(output_folder)
    
def file_reader(file_path):
    dfs = pd.read_csv(file_path)
    return dfs

def concat_meandwell(input_folder):
    filenames = glob.glob(input_folder + "/*.csv")

    mean_concat = []
    for filename in filenames:
        mean_concat.append(pd.read_csv(filename))
    test = pd.concat(mean_concat, ignore_index=True)
    test_dfs = pd.DataFrame(test)
    return test_dfs

hist_data = file_reader(f'{input_folder}/Raw_FRET_histogram_data.csv')
transition_frequency = file_reader(f'{input_folder}/Transition_frequency.csv')
mean_dwell = concat_meandwell(f'{input_folder}/Mean_dwell')

def float_generator(data_frame, treatment):
    transition_frequency_arrow = data_frame[data_frame["sample"]== treatment]
    normalised_number_lowtohigh = float(np.array(transition_frequency_arrow["< 0.5 to > 0.5"])/1000)
    normalised_number_hightolow = float(np.array(transition_frequency_arrow["> 0.5 to < 0.5"])/1000)
    normalised_number_hightohigh = float(np.array(transition_frequency_arrow["> 0.5 to > 0.5"])/1000)
    normalised_number_lowtolow = float(np.array(transition_frequency_arrow["< 0.5 to < 0.5"])/1000)
    arrow_list = [normalised_number_lowtohigh,normalised_number_hightolow,normalised_number_hightohigh,normalised_number_lowtolow]
    return arrow_list


def heatmap_prep(histogram_data, treatment):
    subset_data = histogram_data[histogram_data["treatment_name"]==treatment]
    total = len(subset_data[(subset_data["FRET"] < FRET_thresh) | (subset_data["FRET"] > FRET_thresh)])
    subset_data_largerthanthresh = len(subset_data[subset_data["FRET"] > FRET_thresh])
    subset_data_lessthanthresh = len(subset_data[subset_data["FRET"] < FRET_thresh])
    time_below = (subset_data_lessthanthresh/total)
    time_above = (subset_data_largerthanthresh/total)
    thresh_dicts = {"< 0.5":[time_below], "> 0.5":[time_above] }
    thresh_dfs = pd.DataFrame(thresh_dicts)
    #thresh_dfs["treatment"] = treatment
    return thresh_dfs

def mean_dwell_prep(mean_dwell_data, treatment):
    subset_mean_dwell = mean_dwell_data[mean_dwell_data["sample"]== treatment]
    meandwell_lowtohigh = float(np.array(subset_mean_dwell["< 0.5 to > 0.5"]))
    meandwell_hightolow = float(np.array(subset_mean_dwell["> 0.5 to < 0.5"]))
    meandwell_hightohigh = float(np.array(subset_mean_dwell["> 0.5 to > 0.5"]))
    meandwell_lowtolow = float(np.array(subset_mean_dwell["< 0.5 to < 0.5"]))
    mean_list = [meandwell_lowtohigh,meandwell_hightolow,meandwell_hightohigh,meandwell_lowtolow]
    return mean_list

def plot_heatmap(df_to_plot, arrow_details, mean):
    fig, ax = plt.subplots(figsize=(6,7))
    plt.rcParams['svg.fonttype'] = 'none'
    sns.set(font_scale = 1.5)
    sns.heatmap(df_to_plot, cmap = "PuBuGn", annot=True, linewidths=5, linecolor= "black", cbar_kws={"label":"fraction of total time",'use_gridspec': False, 'location': 'top'}, vmin = 0, vmax = 1)
    ax.arrow(.8,.225, .3,0, width=arrow_details[0], color="darkorange", head_width=arrow_details[0]*4, head_length=0.1)
    ax.arrow(1.2,.775, -.3,0, width=arrow_details[1], color="darkorange", head_width=arrow_details[1]*4, head_length=0.1)
    ax.arrow(1.85,.75,0 ,0.1, width=arrow_details[2], color="darkorange", head_width=arrow_details[2]*2, head_length=0.1)
    ax.arrow(1.65,.95,0 ,-0.1, width=arrow_details[2], color="darkorange", head_width=arrow_details[2]*2, head_length=0.1)
    ax.arrow(0.15,.05,0 ,0.1, width=arrow_details[3], color="darkorange", head_width=arrow_details[3]*2, head_length=0.1)
    ax.arrow(0.35,.25,0,-0.1, width=arrow_details[3], color="darkorange", head_width=arrow_details[3]*2, head_length=0.1)
    plt.text(.69, .3, str(f"{float(mean[0]):.3}" + " s"))
    plt.text(.10, .3, str(f"{float(mean[3]):.3}" + " s"))
    plt.text(1.05, .725, str(f"{float(mean[1]):.3}" + " s"))
    plt.text(1.6, .725, str(f"{float(mean[2]):.3}" + " s"))
    return fig


treatment_list = mean_dwell["sample"].to_list()  ### add .unique() if you were to use the hist_data dataframe

for data_name in treatment_list:
    arrow_info = float_generator(transition_frequency, data_name)
    heatmap_info = heatmap_prep(hist_data, data_name)
    mean = mean_dwell_prep(mean_dwell, data_name)
    plot_heatmap(heatmap_info, arrow_info, mean).savefig(f"{output_folder}/Heatmap_{data_name}.svg", dpi = 600)


