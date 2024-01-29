from scipy.stats import norm
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import pandas as pd
import numpy as np
import glob
import os
from Utilities.Data_analysis import filter_TDP, file_reader, count_filt_mol

output_folder = 'Experiment_1-description/python_results'
plot_export = f'{output_folder}/synchronised_transitions/'
if not os.path.exists(plot_export):
    os.makedirs(plot_export)
compiled_df_HMM = pd.read_csv(f'{output_folder}/compiled_df_HMM.csv')

FRET_value = 0.5 #### the proportion of molecules that travel below this threshold will be counted
exposure = 0.2 # in seconds
order = ['treatment1', 'treatment2']
list_to_drop = ['']
frames_to_plot = 50
FRET_before = 0.5
FRET_after = 0.5

compiled_df_HMM_dropped = compiled_df_HMM.dropna()
###########
########### Calculates the dwell time for each e_pred_global state based on the HMM fits and then appends the dwell duration to the cleaned FRET histogram data.
########### The dwell state (i.e., state 1, state 2, etc) is also appended to the dataset. This code also removes molecules that have only a single dwell state (i.e., does not exhibit a transition)
###########

def calculate_dwells(dfs):
    """Calculate the dwell duration and number for each molecule and then appends it to the cleaned histogram dataset

    Args:
        dfs (df): dataframe containing the cleaned FRET histogram data. This dataset is produced in the 1A-plot_histogram script

    Returns:
        df: returns the same dataframe with addition columns containing the dwell state and the duration of that state
    """
    compiled = []
    for molecule, df in dfs.groupby(['cumulative_molecule']):
        frame_length = len(df[df['cumulative_molecule']== molecule])
        df['dwell_steady_state'] = df['e_pred_global'].ne(df['e_pred_global'].shift()).cumsum()
        test_df = pd.DataFrame(df.groupby([df['e_pred_global'].ne(df['e_pred_global'].shift()).cumsum(), 'e_pred_global']).size())
        test_df.index.names = ["transition", "e_pred_global"]
        test_df.reset_index(inplace = True)
        test_df.columns = ["transition", "e_pred_global", 'dwell']
        dict_dwell = dict(zip(test_df['transition'], test_df['dwell']))
        df['dwell'] = df['dwell_steady_state'].map(dict_dwell)
        df['frame_length'] = frame_length
        compiled.append(df)
    compiled_df = pd.concat(compiled)
    ####### This next section removes any molecules that have only a single dwell state (i.e., they photobleach
    # before any transition occcurs)
    filtered = [df2 for (molecule, treatment), df2 in compiled_df.groupby(['cumulative_molecule', 'treatment']) if df2['dwell_steady_state'].nunique() > 1]
    return pd.concat(filtered)


##################
##################  Identifies the transition point (the point at which the e_pred_global changes), and then assigns that transition a dwell time and the FRET states before and 
##################  after a transition.
##################

def generate_transitions(dfs):
    """identifies the time at which a transition occurs and provides the FRET state before (FRET_before) and after (FRET_after) a transition occurs.

    Args:
        dfs (df): dataframe containing the cleaned FRET histogram data with the dwell time of each state for each molecule. Generated using the calculate dwells function.

    Returns:
        df: dataframe that contains extra columns, which include the transition point (i.e., the point at which the e_pred_global changes, is either True or False), transition dwell (the duration of FRET_before prior to a True transition) and the FRET_before or FRET_after a transition
    """
    compiled_transition = []
    for molecule, df in dfs.groupby(['cumulative_molecule']):
        dwell_steady_state_list = list(df['dwell_steady_state'].unique())
        df['transition_point'] = df['e_pred_global'].ne(df['e_pred_global'].shift())
        df['transition_point'].iloc[0] = False
        df['column_for_dict'] = df['dwell_steady_state']-1
        steady_dwell = df[['dwell_steady_state', 'dwell']]
        dwell_dict = test_dict = dict(zip(steady_dwell['dwell_steady_state'], steady_dwell['dwell']))
        df['transition_dwell'] = df['column_for_dict'].map(dwell_dict)
        steadyFRET = df[['dwell_steady_state', 'e_pred_global']]
        test_dict = dict(zip(steadyFRET['dwell_steady_state'], steadyFRET['e_pred_global']))
        df['FRET_before'] = df['column_for_dict'].map(test_dict)
        # df.dropna('index', inplace = True)
        df['FRET_after'] = df['e_pred_global']
        df.drop('column_for_dict', axis = 1, inplace = True)
        compiled_transition.append(df)
    return pd.concat(compiled_transition)



#########
######### Filters the dataset prior to plotting. Will filter based on the type of transition (low-to-high or vice-versa) and for how long the dwell state exists prior to a transition.
######### The next function will then plot all the transitions that meet the filtered criteria.
#########

def filt_df_to_plot(df, FRET_before, FRET_after, transition_type = 'low_to_high', min_dwell_before = 0):
    """will filter the dataframe according to the transition of interest and the dwell time of the FRET state prior to that transition. Returns a list of indexes that meet the transition criteria

    Args:
        df (df): dataframe containing the cleaned FRET data with transition information
        FRET_before (float): FRET state prior to transition, used to filter data
        FRET_after (floar): FRET state after transtion, used to filter data
        transition_type (str, optional): determines what kind of transitions you want to look into (e.g., low-to-high transitions where low is below FRET_before and high is above FRET_after). Defaults to 'low_to_high'.
        min_dwell_before (int, optional): variable that defines for how long a FRET state must have existed before the transition. Defaults to 0.

    Returns:
        list: returns a list of index values where the above transition criteria is true. This list is then used to identify transition points within the cleaned histogram data and plot.
    """
    transitions_to_plot = df[df['transition_point'] == True]
    if transition_type == 'low_to_high':
        index_to_plot = transitions_to_plot[((transitions_to_plot['FRET_before'] < FRET_before) & (transitions_to_plot['transition_dwell'] > min_dwell_before)) & (transitions_to_plot['FRET_after'] > FRET_after)].index
    elif transition_type == 'high_to_low':
        index_to_plot = transitions_to_plot[((transitions_to_plot['FRET_before'] > FRET_before) & (transitions_to_plot['transition_dwell'] > min_dwell_before)) & (transitions_to_plot['FRET_after'] < FRET_after)].index
    return index_to_plot


def plot_synchronised_transition(dfs, index_to_plot, exposure_seconds, list_to_drop, order, frame_from_trans = 80, label = ''):
    """plots the FRET values either side of a transition type of interest

    Args:
        dfs (df): dataframe containing the raw FRET values (generated after the calculate dwells function)
        index_to_plot (list): list of index values that met the criteria defined in 'filt_df_to_plot' that will the be mapped to df for plotting of raw FRET values
        exposure_seconds (float): exposure in seconds used to convert frames to a unit of time
        list_to_drop (list): list containing the name of all treatments to be dropped
        frame_from_trans (int, optional): should be the same as 'min_dwell_before' variable in the 'filt_df_to_plot' function, basically sets the xlims. Defaults to 80.
    """
    combined_mini = []
    for df in index_to_plot:
        print(df)
        lower = df - frame_from_trans
        upper = df + (frame_from_trans+1)
        mini_df = dfs.iloc[lower:upper].reset_index()
        mini_df['time_from_trans'] = np.arange(-(frame_from_trans/(1/exposure_seconds)), (frame_from_trans/(1/exposure_seconds))+exposure_seconds, exposure_seconds)
        combined_mini.append(mini_df)
    combined_mini = pd.concat(combined_mini)
    filt_data = combined_mini[~combined_mini['treatment'].isin(list_to_drop)]
    fig, axes = plt.subplots()
    sns.set(style = 'ticks')
    sns.lineplot(data = filt_data, x = 'time_from_trans', y = 'FRET Cy3 to AF647', hue = 'treatment', palette = 'BuPu', hue_order = order)
    plt.xlabel('Time (s)')
    plt.legend(title = '',loc = 'best')
    fig.savefig(f'{plot_export}/synchronised_release{"_"+label}.svg', dpi = 600)
    plt.show()

calculated_transitions = []
for treatment, df in compiled_df_HMM_dropped.groupby('treatment'):
    dwell_df = calculate_dwells(df)
    transition_df = generate_transitions(dwell_df)
    calculated_transitions.append(transition_df)
calculated_transitions_df = pd.concat(calculated_transitions)
calculated_transitions_df.reset_index(inplace = True)
calculated_transitions_df.drop('index', axis = 1, inplace = True)

font = {'weight' : 'normal', 'size'   : 12 }
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.family'] = "sans-serif"
plt.rc('font', **font)
plt.rcParams['svg.fonttype'] = 'none'

dnak_stable_release = filt_df_to_plot(calculated_transitions_df, FRET_before, FRET_after,'low_to_high', frames_to_plot)
plot_synchronised_transition(calculated_transitions_df, dnak_stable_release, exposure, list_to_drop, order, frames_to_plot, 'release')



def plot_synchronised_fluorescence(dfs, index_to_plot, exposure_seconds, list_to_drop, frame_from_trans = 80, label = ''):
    """plots the FRET values and total fluorescence of all dyes following excitation at 488 nm either side of a transition type of interest

    Args:
        dfs (df): dataframe containing the raw FRET values (generated after the calculate dwells function)
        index_to_plot (list): list of index values that met the criteria defined in 'filt_df_to_plot' that will the be mapped to df for plotting of raw FRET values
        exposure_seconds (float): exposure in seconds used to convert frames to a unit of time
        list_to_drop (list): list containing the name of all treatments to be dropped
        frame_from_trans (int, optional): should be the same as 'min_dwell_before' variable in the 'filt_df_to_plot' function, basically sets the xlims. Defaults to 80.
    """
    combined_mini = []
    for df in index_to_plot:
        print(df)
        lower = df - frame_from_trans
        upper = df + (frame_from_trans+1)
        mini_df = dfs.iloc[lower:upper].reset_index()
        mini_df['time_from_trans'] = np.arange(-(frame_from_trans/(1/exposure_seconds)), (frame_from_trans/(1/exposure_seconds))+exposure_seconds, exposure_seconds)
        combined_mini.append(mini_df)
    combined_mini = pd.concat(combined_mini)
    filt_data = combined_mini[~combined_mini['treatment'].isin(list_to_drop)]
    for treatment, df in filt_data.groupby('treatment'):
        fig, axes = plt.subplots()
        sns.set(style = 'ticks')
        ax2 = axes.twinx()
        sns.lineplot(data = df, x = 'time_from_trans', y = 'FRET Cy3 to AF647', color = 'black', ax = axes)
        sns.lineplot(data = df, x = 'time_from_trans', y = 'normalised_summed_fluorescence', color = 'skyblue', ax = ax2)
        axes.set_xlabel('Time (s)')
        ax2.set_ylabel('Normalised total fluorescence (a.u.)')  
        axes.set_ylabel('FRET')  
        axes.legend(['FRET'], loc = 'upper left')
        ax2.legend(['Fluorescence'])
        fig.savefig(f'{plot_export}/{treatment}_synchronised_release_fluorescence{"_"+label}.svg', dpi = 600)
        plt.show()

bottle = []
for (treatment, molecule), df in calculated_transitions_df.groupby(['treatment', 'cumulative_molecule']):
    df['normalised_summed_fluorescence'] = df['probe_summed_fluorescence']/df['probe_summed_fluorescence'].max()
    bottle.append(df)
calculated_transitions_df = pd.concat(bottle)


dnak_stable_release = filt_df_to_plot(calculated_transitions_df, FRET_before, FRET_after,'low_to_high', frames_to_plot)
dnak_stable_binding = filt_df_to_plot(calculated_transitions_df, FRET_before, FRET_after, 'high_to_low', frames_to_plot)

plot_synchronised_fluorescence(calculated_transitions_df, dnak_stable_release, exposure, list_to_drop, frames_to_plot, 'release')
plot_synchronised_fluorescence(calculated_transitions_df, dnak_stable_binding, exposure, list_to_drop, frames_to_plot, 'binding')
