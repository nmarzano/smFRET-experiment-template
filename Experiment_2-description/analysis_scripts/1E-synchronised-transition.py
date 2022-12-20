from scipy.stats import norm
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import pandas as pd
import numpy as np
import glob
import os
from Utilities.Data_analysis import filter_TDP, file_reader, count_filt_mol

output_folder = 'Experiment/python_results'
compiled_data = pd.read_csv(f'{output_folder}/Cleaned_FRET_histogram_data.csv')
FRET_value = 0.5 #### the proportion of molecules that travel below this threshold will be counted
exposure = 0.2 # in seconds

def calculate_dwell2(dfs):
    compiled = []
    for molecule, df in dfs.groupby(['molecule_number']):
        frame_length = len(df[df['molecule_number']== molecule])
        df['dwell_steady_state'] = df['idealized FRET'].ne(df['idealized FRET'].shift()).cumsum()
        test_df = pd.DataFrame(df.groupby([df['idealized FRET'].ne(df['idealized FRET'].shift()).cumsum(), 'idealized FRET']).size())
        test_df.index.names = ["transition", "idealized FRET"]
        test_df.reset_index(inplace = True)
        test_df.columns = ["transition", "idealized FRET", 'dwell']
        dict_dwell = dict(zip(test_df['transition'], test_df['dwell']))
        df['dwell'] = df['dwell_steady_state'].map(dict_dwell)
        df['frame_length'] = frame_length
        compiled.append(df)
    compiled_df = pd.concat(compiled)
    ####### This next section removes any molecules that have only a single dwell state (i.e., they photobleach
    # before any transition occcurs)
    filtered = []
    for (molecule, treatment), df2 in compiled_df.groupby(['molecule_number', 'treatment_name']):
        if df2['dwell_steady_state'].nunique() > 1:
            filtered.append(df2)
    return pd.concat(filtered)

def generate_transitions2(dfs):
    compiled_transition = []
    for molecule, df in dfs.groupby(['molecule_number']):
        dwell_steady_state_list = list(df['dwell_steady_state'].unique())
        df['transition_point'] = df['idealized FRET'].ne(df['idealized FRET'].shift())
        df['transition_point'].iloc[0] = False
        df['column_for_dict'] = df['dwell_steady_state']-1
        steady_dwell = df[['dwell_steady_state', 'dwell']]
        dwell_dict = test_dict = dict(zip(steady_dwell['dwell_steady_state'], steady_dwell['dwell']))
        df['transition_dwell'] = df['column_for_dict'].map(dwell_dict)
        steadyFRET = df[['dwell_steady_state', 'idealized FRET']]
        test_dict = dict(zip(steadyFRET['dwell_steady_state'], steadyFRET['idealized FRET']))
        df['FRET_before'] = df['column_for_dict'].map(test_dict)
        # df.dropna('index', inplace = True)
        df['FRET_after'] = df['idealized FRET']
        df.drop('column_for_dict', axis = 1, inplace = True)
        compiled_transition.append(df)
    return pd.concat(compiled_transition)


poop = []
for treatment, df in compiled_data.groupby('treatment_name'):
    data = calculate_dwell2(df)
    data2 = generate_transitions2(data)
    poop.append(data2)
pooped = pd.concat(poop)
pooped.reset_index(inplace = True)
pooped.drop('index', axis = 1, inplace = True)


def filt_df_to_plot(df, FRET_before, FRET_after, transition_type = 'low_to_high', min_dwell_before = 0):
    transitions_to_plot = df[df['transition_point'] == True]
    if transition_type == 'low_to_high':
        index_to_plot = transitions_to_plot[((transitions_to_plot['FRET_before'] < FRET_before) & (transitions_to_plot['transition_dwell'] > min_dwell_before)) & (transitions_to_plot['FRET_after'] > FRET_after)].index
    elif transition_type == 'high_to_low':
        index_to_plot = transitions_to_plot[((transitions_to_plot['FRET_before'] > FRET_before) & (transitions_to_plot['transition_dwell'] > min_dwell_before)) & (transitions_to_plot['FRET_after'] < FRET_after)].index
    return index_to_plot


def plot_synchronised_transition(dfs, index_to_plot, exposure_seconds, list_to_drop, frame_from_trans = 80):
    combined_mini = []
    for df in index_to_plot:
        print(df)
        lower = df - frame_from_trans
        upper = df + (frame_from_trans+1)
        mini_df = dfs.iloc[lower:upper].reset_index()
        mini_df['time_from_trans'] = np.arange(-(frame_from_trans/(1/exposure_seconds)), (frame_from_trans/(1/exposure_seconds))+exposure_seconds, exposure_seconds)
        combined_mini.append(mini_df)
    combined_mini = pd.concat(combined_mini)
    filt_data = combined_mini[~combined_mini['treatment_name'].isin(list_to_drop)]
    fig, axes = plt.subplots()
    sns.set(style = 'ticks')
    sns.lineplot(data = filt_data, x = 'time_from_trans', y = 'FRET', hue = 'treatment_name')
    plt.xlabel('Time (s)')
    plt.legend(loc = 'best')
    plt.show()


list_to_drop = ['PPR_alone', 'RNA05']
dnak_stable_release = filt_df_to_plot(pooped, 0.5, 0.5,'low_to_high', 30)
plot_synchronised_transition(pooped, dnak_stable_release, exposure, list_to_drop, 30)

dnak_stable_binding = filt_df_to_plot(pooped, 0.5, 0.5, 'high_to_low', 0)
plot_synchronised_transition(pooped, dnak_stable_binding, exposure, list_to_drop, 100)

