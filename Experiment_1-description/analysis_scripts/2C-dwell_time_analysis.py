from scipy.stats import norm
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import pandas as pd
import numpy as np
import glob

output_folder = 'Experiment_X-description/python_results'
filename = f'{output_folder}/TDP_cleaned.csv'


FRET_thresh = 0.5 #### FRET value at which to filter data above or below. IF CHANGED, WILL NEED TO CHANGE ALL 0.5 VALUES (E.G. BELOW IN HEADERS) TO THE NEW VALUE
fps = 5  ### frames per second
thresh = 2 ### should be 10x expsoure if using NL515 smoothing on MASH FRET
headers = [f"< {FRET_thresh} to < {FRET_thresh}", f"< {FRET_thresh} to > {FRET_thresh}", f"> {FRET_thresh} to > {FRET_thresh}", f"> {FRET_thresh} to < {FRET_thresh}"]
headers_withsum =  [f"< {FRET_thresh} to < {FRET_thresh}", f"< {FRET_thresh} to > {FRET_thresh}", f"> {FRET_thresh} to > {FRET_thresh}", f"> {FRET_thresh} to < {FRET_thresh}", "> 0.5 to < 0.5", "sum", "sample"]
TDP_data = pd.read_csv(filename, header = "infer")

from Utilities.Data_analysis import cleanup_dwell, filter_dwell, transition_frequency, calculate_mean, fret_before_trans


for treatment_name, df in TDP_data.groupby("treatment_name"):
    initial_data = df[df["treatment_name"] == treatment_name]
    cleaned_data = cleanup_dwell(initial_data, fps, thresh) ##### to keep the first dwell state, simply change code to "cleanup_dwell(initial_data, "keep")
    filtered_data = filter_dwell(cleaned_data, FRET_thresh, headers)
    filtered_data.to_csv(f"{output_folder}/Dwell_times/Filtered_dwelltime_{treatment_name}.csv", index = False)
    mean_dwell = calculate_mean(filtered_data, treatment_name)
    mean_dwell.to_csv(f"{output_folder}/Mean_dwell/Filtered_meandwell_{treatment_name}.csv", index = False)
    dwell_frequency = transition_frequency(filtered_data)
    dwell_frequency["sample"] = treatment_name
    dwell_frequency.to_csv(f"{output_folder}/Dwell_frequency/Filtered_dwellfrequency_{treatment_name}.csv", index = False, header = None)


FRET_value_before_transition = fret_before_trans(TDP_data, 0.7, fps, FRET_thresh)

def plot_fret_before(df, to_drop = 'none'):
    if to_drop == 'none':
        plot1 = plt.figure(figsize = (12, 6))
        sns.set(style = "darkgrid", font_scale = 1.5)
        sns.violinplot(data = df, x = 'treatment_name', y = 'FRET before transition')
        sns.stripplot(data = df, x = 'treatment_name', y = 'FRET before transition', color='black', alpha = 0.5)
    else:
        dropped = df[~df['treatment_name'].isin(to_drop)].dropna()
        plot1 = plt.figure(figsize = (12, 6))
        sns.set(style = "darkgrid", font_scale = 1.5)
        sns.violinplot(data = dropped, x = 'treatment_name', y = 'FRET before transition')
        sns.stripplot(data = dropped, x = 'treatment_name', y = 'FRET before transition', color='black', alpha = 0.5)
    plt.rcParams['svg.fonttype'] = 'none'
    plt.xlabel('Treatment')
    plt.ylabel('FBefore')
    plot1.savefig(f'{output_folder}/FRET_before_trans.svg', dpi = 600)
    return plot1

plot_fret_before(FRET_value_before_transition)


