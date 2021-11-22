from scipy.stats import norm
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import pandas as pd
import numpy as np
import glob
import functools
import os

output_folder = 'Experiment_X-description/python_results'
plot_folder = f'{output_folder}/dwell_analysis_figs'
if not os.path.exists(plot_folder):
    os.makedirs(plot_folder)
filename = f'{output_folder}/TDP_cleaned.csv'

order = ['Native', 'Spontaneous', 'KJ', 'low_GrpE', 'high_GrpE']
palette = 'mako'
FRET_thresh = 0.5 #### FRET value at which to filter data above or below. IF CHANGED, WILL NEED TO CHANGE ALL 0.5 VALUES (E.G. BELOW IN HEADERS) TO THE NEW VALUE
fps = 5  ### frames per second
thresh = 2 ### should be 10x expsoure if using NL515 smoothing on MASH FRET
headers = [f"< {FRET_thresh} to < {FRET_thresh}", f"< {FRET_thresh} to > {FRET_thresh}", f"> {FRET_thresh} to > {FRET_thresh}", f"> {FRET_thresh} to < {FRET_thresh}"]
headers_withsum =  [f"< {FRET_thresh} to < {FRET_thresh}", f"< {FRET_thresh} to > {FRET_thresh}", f"> {FRET_thresh} to > {FRET_thresh}", f"> {FRET_thresh} to < {FRET_thresh}", "> 0.5 to < 0.5", "sum", "sample"]
TDP_data = pd.read_csv(filename, header = "infer")

from Utilities.Data_analysis import cleanup_dwell, filter_dwell, transition_frequency, calculate_mean, fret_state_trans


for treatment_name, df in TDP_data.groupby("treatment_name"):
    initial_data = df[df["treatment_name"] == treatment_name]
    cleaned_data = cleanup_dwell(initial_data, fps, thresh, 'keep') ##### to keep the first dwell state, simply change code to "cleanup_dwell(initial_data, "keep")
    filtered_data = filter_dwell(cleaned_data, FRET_thresh, headers)
    filtered_data.to_csv(f"{output_folder}/Dwell_times/Filtered_dwelltime_{treatment_name}.csv", index = False)
    mean_dwell = calculate_mean(filtered_data, treatment_name)
    mean_dwell.to_csv(f"{output_folder}/Mean_dwell/Filtered_meandwell_{treatment_name}.csv", index = False)
    dwell_frequency = transition_frequency(filtered_data)
    dwell_frequency["sample"] = treatment_name
    dwell_frequency.to_csv(f"{output_folder}/Dwell_frequency/Filtered_dwellfrequency_{treatment_name}.csv", index = False, header = None)


Transition_threshold = 0.5

def plot_fret_trans(df, FRET_state = 'after', to_drop = 'none', threshold = Transition_threshold):
    if to_drop == 'none':
        if FRET_state == 'after':
            plot1 = plt.figure(figsize = (12, 6))
            sns.set(style = "darkgrid", font_scale = 1.5)
            sns.violinplot(data = df, x = 'treatment_name', y = 'FRET_after', palette = palette, order = order)
            sns.stripplot(data = df, x = 'treatment_name', y = 'FRET_after', color='black', alpha = 0.25, order = order)
            plt.ylabel(f'FRET state after transition from < {threshold}')
        elif FRET_state == 'before':
            plot1 = plt.figure(figsize = (12, 6))
            sns.set(style = "darkgrid", font_scale = 1.5)
            sns.violinplot(data = df, x = 'treatment_name', y = 'FRET_before', palette = palette, order = order)
            sns.stripplot(data = df, x = 'treatment_name', y = 'FRET_before', color='black', alpha = 0.25, order = order)
            plt.ylabel(f'FRET state before transition to < {threshold}')
    else:
        dropped = df[~df['treatment_name'].isin(to_drop)].dropna()
        plot1 = plt.figure(figsize = (12, 6))
        sns.set(style = "darkgrid", font_scale = 1.5)
        sns.violinplot(data = dropped, x = 'treatment_name', y = 'FRET_before')
        sns.stripplot(data = dropped, x = 'treatment_name', y = 'FRET_before', color='black', alpha = 0.25)
    plt.rcParams['svg.fonttype'] = 'none'
    plt.xlabel('Treatment')
    plt.ylim(-0.1, 1.2)
    plt.xticks(rotation=45)
    plot1.savefig(f'{plot_folder}/FRET_{FRET_state}_trans_{Transition_threshold}.svg', dpi = 600)
    plt.show()


FRET_value_after_transition = fret_state_trans(TDP_data, Transition_threshold, fps, FRET_thresh, 'after')
plot_fret_trans(FRET_value_after_transition, 'after')


FRET_value_before_transition = fret_state_trans(TDP_data, Transition_threshold, fps, FRET_thresh, 'before')
plot_fret_trans(FRET_value_before_transition, 'before')



def count_chaperone_events(dfs, thresh, fps_clean, thresh_clean):
    cleaned_df = []
    for treatment_name, df in dfs.groupby("treatment_name"):
        initial_data = df[df["treatment_name"] == treatment_name]    
        cleaned = cleanup_dwell(initial_data, fps_clean, thresh_clean)
        cleaned_df.append(cleaned)
    cleaned_concat = pd.concat(cleaned_df)
    cleaned_concat['Total Molecule Lifetime (min)'] = (cleaned_concat['number_of_frames']/5)/60
    filt = []
    for treatment_name, df in cleaned_concat.groupby("treatment_name"):
        treatment = treatment_name
        chaperone_on = df[(df['FRET_after'] <= thresh) & (df['FRET_before'] >= thresh)].groupby('Molecule').count()['FRET_after'].reset_index()
        chaperone_off = df[(df['FRET_after'] >= thresh) & (df['FRET_before'] <= thresh)].groupby('Molecule').count()['FRET_before'].reset_index()
        time = df.groupby('Molecule').mean()['Total Molecule Lifetime (min)'].reset_index()
        # merged_test = chaperone_on.merge(chaperone_off, how = 'outer').fillna(0)
        merged_test = functools.reduce(lambda left, right: pd.merge(left, right, on='Molecule', how='outer'), [chaperone_on, chaperone_off, time]) ### Really usefull code for merging multuple dfs
        merged_test['treatment'] = treatment
        filt.append(merged_test)
    count_data = pd.concat(filt)
    test = pd.DataFrame(count_data)
    test.dropna(subset = ['FRET_after', 'FRET_before'], how = 'all', inplace = True)
    test.fillna(0, inplace = True)
    return test
 
org_chap_events = count_chaperone_events(dfs = TDP_data, thresh = 0.3, fps_clean = fps, thresh_clean = FRET_thresh)
org_chap_events['FRET_after_normalised'] = org_chap_events['FRET_after']/org_chap_events['Total Molecule Lifetime (min)']
org_chap_events['FRET_before_normalised'] = org_chap_events['FRET_before']/org_chap_events['Total Molecule Lifetime (min)']

def plot_binding_release(df, chaperone = 'binding'):
    if chaperone == 'binding':
        plot1 = plt.figure(figsize = (12, 6))
        sns.violinplot(data = df, y = 'FRET_after_normalised', x = 'treatment', cut = 0, order = order, palette = palette)
        sns.stripplot(data = df, y = 'FRET_after_normalised', x = 'treatment', color='black', alpha = 0.25, order = order)
        plt.ylabel('# of chaperone binding events/min/molecule')
        plt.xticks(rotation = 45)
        plot1.savefig(f'{plot_folder}/chaperone_binding_rate_per_molecule.svg', dpi = 600)
        plt.show()
    if chaperone == 'release':
        plot2 = plt.figure(figsize = (12, 6))
        sns.violinplot(data = df, y = 'FRET_before_normalised', x = 'treatment', cut = 0, order = order, palette = palette)
        sns.stripplot(data = df, y = 'FRET_before_normalised', x = 'treatment', color='black', alpha = 0.25, order = order)
        plt.ylabel('# of chaperone release events/min/molecule')
        plt.xticks(rotation = 45)
        plot2.savefig(f'{plot_folder}/chaperone_release_rate_per_molecule.svg', dpi = 600)
        plt.show() 
    if chaperone == 'binding_events':
        plot3 = plt.figure(figsize = (12, 6))
        sns.violinplot(data = df, y = 'FRET_after', x = 'treatment', cut = 0, order = order, palette = palette)
        sns.stripplot(data = df, y = 'FRET_after', x = 'treatment', color='black', alpha = 0.25, order = order)
        plt.xticks(rotation = 45)
        plt.ylabel('# of chaperone binding events/molecule')
        plot3.savefig(f'{plot_folder}/chaperone_binding_events_per_molecule.svg', dpi = 600)
        plt.show()

plot_binding_release(org_chap_events, 'binding')  ######## options are 'binding', 'release', 'binding_events'

