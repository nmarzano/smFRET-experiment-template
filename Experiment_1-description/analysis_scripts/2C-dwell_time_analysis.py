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

###############
###############   Plot FRET states before or after a transition to a defined FRET state
###############
Transition_threshold = 0.5

def plot_fret_trans(df, FRET_state = 'after', to_drop = 'none', threshold = Transition_threshold, palette = 'mako'):
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

###############
###############  Calculate the number of binding or release events (defined when FRET crosses a FRET threshold) for each
###############  molecule then normalise to the lifetime of that molecule to get the rate - also plot

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
 
org_chap_events = count_chaperone_events(dfs = TDP_data, thresh = 0.2, fps_clean = fps, thresh_clean = FRET_thresh)
org_chap_events['FRET_after_normalised'] = org_chap_events['FRET_after']/org_chap_events['Total Molecule Lifetime (min)']
org_chap_events['FRET_before_normalised'] = org_chap_events['FRET_before']/org_chap_events['Total Molecule Lifetime (min)']
org_chap_events['bind_and_release'] = org_chap_events[['FRET_after', 'FRET_before']].min(axis = 1)
org_chap_events['bind_and_release_overtime'] = (org_chap_events['bind_and_release']/org_chap_events['Total Molecule Lifetime (min)'])
org_chap_events['bind_and_release_overtime'] = org_chap_events['bind_and_release_overtime'].replace(0, np.nan)


def plot_binding_release(df, chaperone = 'binding', order = False, palette = 'mako'):
    if chaperone == 'binding':
        ycol = 'FRET_after_normalised'
        ylabel = '# of chaperone binding events/min/molecule'
        title = 'chaperone_binding_rate_per_molecule'
    if chaperone == 'release':
        ycol = 'FRET_before_normalised'
        ylabel = '# of chaperone release events/min/molecule'
        title = 'chaperone_release_rate_per_molecule'
    if chaperone == 'binding_events':
        ycol = 'FRET_after'
        ylabel = '# of chaperone binding events/molecule'
        title = 'chaperone_binding_events_per_molecule'
    if chaperone == 'binding_and_release':
        ycol = 'bind_and_release_overtime'
        ylabel = '# of chaperone binding and release events/molecule/min'
        title = 'chaperone_binding_and_release_events_per_molecule_min'
    plot1 = plt.figure(figsize = (12, 6))
    sns.violinplot(data = df, y = ycol, x = 'treatment', cut = 0, order = order, palette = palette)
    sns.stripplot(data = df, y = ycol, x = 'treatment', color='black', alpha = 0.25, order = order)
    plt.xticks(rotation = 45)
    plt.ylabel(f'{ylabel}')
    plot1.savefig(f'{plot_folder}/{title}.svg', dpi = 600)
    plt.show()

plot_binding_release(org_chap_events, 'binding_and_release', order)

#################
#################   Calculate the mean, sem and N of binding and release events for statistical analysis
#################

mean_chaperone_bind_release_mean = org_chap_events.groupby('treatment')['bind_and_release_overtime'].mean()
mean_chaperone_bind_release_sem = org_chap_events.groupby('treatment')['bind_and_release_overtime'].sem()
mean_chaperone_bind_release_N = org_chap_events.groupby('treatment')['bind_and_release_overtime'].count()
bind_release_col = pd.concat([mean_chaperone_bind_release_mean,mean_chaperone_bind_release_sem, mean_chaperone_bind_release_N], axis = 1)
bind_release_col.columns = ['mean_bind_and_release', 'sem', 'n']
bind_release_col.reset_index(inplace = True)
bind_release_col.to_csv(f"{output_folder}/bind_release_col.csv", index = False)

##################
##################  Caluclate the proportion of transitions that are significant (i.e., have a change in FRET state
##################  above a defined threshold) and also the proportion of molecules that experience these transitions
##################
TDP_data['FRET_trans_difference'] = abs(TDP_data['FRET_before'] - TDP_data['FRET_after'])

def find_large_transitions(dfs, delta_thresh):
    mol_with_large_trans = []
    for treatment, df in dfs.groupby('treatment_name'):
        filt = df[df['FRET_trans_difference'] > delta_thresh]
        filt_count_mol = pd.DataFrame(filt[filt['FRET_trans_difference'] > delta_thresh].agg({"Molecule": "nunique"})/df.agg({"Molecule": "nunique"})*100)
        filt_count_mol['treatment'] = treatment
        filt_count_mol['proportion_of_mol'] = (filt['Molecule'].count()/df['Molecule'].count())*100
        mol_with_large_trans.append(filt_count_mol)
    dfs = pd.concat(mol_with_large_trans)
    dfs = pd.DataFrame(dfs).reset_index()
    dfs.drop(columns = 'index', inplace =  True)
    dfs.columns = ['proportion_mol_large_transition', 'treatment', 'proportion_of_large_transitions']
    return dfs

def plot_large_transitions(df, type = 'transition_prob', palette = 'mako'):
    if type == 'transition_prob': 
        ycol = 'proportion_of_large_transitions'
    if type == 'proportion_of_mol': 
        ycol = 'proportion_mol_large_transition'
    plot = plt.figure(figsize = (6, 3))
    sns.barplot(data = df, 
        y = ycol, 
        x = 'treatment',
        order = order, 
        palette = palette)
    plt.xticks(rotation = 45)
    plot.savefig(f'{plot_folder}/{ycol}.svg', dpi = 600)
    plt.show()

large_transitions_to_plot = find_large_transitions(TDP_data, 0.5)
plot_large_transitions(large_transitions_to_plot)      ##### 'proportion_of_mol'

