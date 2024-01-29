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
compiled_data = pd.read_csv(f'{output_folder}/Cleaned_FRET_histogram_data.csv')
FRET_value = 0.5 #### the proportion of molecules that travel below this threshold will be counted
exposure = 0.2 # imaging exposure in seconds
order = ['treatment1', 'treatment2', 'treatment3', 'treatment4']
list_to_drop = [''] ##### add treatments here that you want to exclude from analysis
frames_to_plot = 50  ###### This value is used for two things. (1) This is the number of frames either side of the transition point that will be plotted. (2) This is the minimum
########################### number of frames that a FRET state must exist prior to a transition in order for it to be plotted (this will generally exclude noisy molecules).
FRET_before = 0.5 #### FRET state before a transition
FRET_after = 0.5 ##### FRET state after the transition

###########
########### Calculates the dwell time for each idealized FRET state based on the HMM fits and then appends the dwell duration to the cleaned FRET histogram data.
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
    filtered = [df2 for (molecule, treatment), df2 in compiled_df.groupby(['molecule_number', 'treatment_name']) if df2['dwell_steady_state'].nunique() > 1]
    return pd.concat(filtered)


##################
##################  Identifies the transition point (the point at which the idealized FRET changes), and then assigns that transition a dwell time and the FRET states before and 
##################  after a transition.
##################

def generate_transitions(dfs):
    """identifies the time at which a transition occurs and provides the FRET state before (FRET_before) and after (FRET_after) a transition occurs.

    Args:
        dfs (df): dataframe containing the cleaned FRET histogram data with the dwell time of each state for each molecule. Generated using the calculate dwells function.

    Returns:
        df: dataframe that contains extra columns, which include the transition point (i.e., the point at which the idealized FRET changes, is either True or False), transition dwell (the duration of FRET_before prior to a True transition) and the FRET_before or FRET_after a transition
    """
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
    filt_data = combined_mini[~combined_mini['treatment_name'].isin(list_to_drop)]
    fig, axes = plt.subplots()
    sns.set(style = 'ticks')
    sns.lineplot(data = filt_data, x = 'time_from_trans', y = 'FRET', hue = 'treatment_name', palette = 'BuPu', hue_order = order)
    plt.xlabel('Time (s)')
    plt.legend(title = '',loc = 'best')
    fig.savefig(f'{plot_export}/synchronised_release{"_"+label}.svg', dpi = 600)
    plt.show()


def determine_first_transition_in_sequence(dfs):
    combined_data = []
    for molecule, df in dfs.groupby('molecule_number'):
        trans_list = df['FRET_increase'].to_list()
        df['output_column'] = [
        1 if all([(x == 0), (trans_list[x] == True), (trans_list[x+1] == True)]) else 
        (1 if all([(trans_list[x-1] == False),(trans_list[x] == True), (trans_list[x+1] == True)]) else 0)
        for x in range(len(trans_list)-1)]+[0]
        combined_data.append(df)
    return pd.concat(combined_data)

def concat_trans_proportion(dfs, raw_df, FRET_before, FRET_after):
    #### create dataframes filtered to have only consecutive or non-consecutive transitions
    consecutive_trans = dfs[dfs['output_column']==True]
    nonconsecutive_trans = dfs[dfs['output_column']==False]

    ##### determine the total number of transitions and the number of transitions that meet a criteria
    transitions_only = raw_df[raw_df['transition_point'] == True]
    total_trans = transitions_only.groupby('treatment_name')['transition_point'].sum()
    transitions_above_thresh = transitions_only[(transitions_only['FRET_before'] < FRET_before) & (transitions_only['FRET_after'] > FRET_after)]
    trans_above_thresh = transitions_above_thresh.groupby('treatment_name')['transition_point'].sum()
    percent_trans_meet_criteria = (trans_above_thresh/total_trans)*100
    percent_trans_meet_criteria_df = pd.DataFrame(percent_trans_meet_criteria).reset_index()
    percent_trans_meet_criteria_df.columns = ['treatment', '% trans DnaK release']

    #### determine the total number of consecutive transitions that meet a criteria
    number_consecutive_event = consecutive_trans.groupby('treatment_name')['output_column'].sum()
    consecutive_event_above_thresh = consecutive_trans[(consecutive_trans['FRET_before'] < FRET_before) & (consecutive_trans['FRET_after'] > FRET_after)]
    number_consecutive_event_meet_criteria = consecutive_event_above_thresh.groupby('treatment_name')['output_column'].sum()
    percent_of_consecutive_event_from_DnaK = (number_consecutive_event_meet_criteria/number_consecutive_event)*100
    percent_of_consecutive_event_from_DnaK = pd.DataFrame(percent_of_consecutive_event_from_DnaK).reset_index()

    percent_of_DnaK_release_events_that_are_consecutive = pd.DataFrame((number_consecutive_event_meet_criteria/trans_above_thresh)*100).reset_index()
    percent_of_DnaK_release_events_that_are_consecutive.columns = ['treatment', 'proportion_consecutive_from_DnaK']

    ### merge columns
    percent_trans_meet_criteria_df['% DnaK release are consecutive'] = percent_of_DnaK_release_events_that_are_consecutive['proportion_consecutive_from_DnaK']
    percent_trans_meet_criteria_df['% consecutive events are DnaK release'] = percent_of_consecutive_event_from_DnaK['output_column']
    return consecutive_trans, nonconsecutive_trans, percent_trans_meet_criteria_df


calculated_transitions = []
for treatment, df in compiled_data.groupby('treatment_name'):
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

dnak_stable_binding = filt_df_to_plot(calculated_transitions_df, FRET_before, FRET_after, 'high_to_low', frames_to_plot)
plot_synchronised_transition(calculated_transitions_df, dnak_stable_binding, exposure, list_to_drop, order, frames_to_plot, 'binding')


############
############ Code past this point is to investigate the difference between consecutive and non-consecutive transitions. At the moment, this code only looks at 
############ consecutive or non-consecutive low-to-high transitions. This can be easily changed by modifing the 'low_to_high' variables below to 'high_to_low'.
############

col = []
for treatment, df in calculated_transitions_df.groupby('treatment_name'):
    transition_data = df[df['transition_point']==True]
    transition_data['FRET_increase'] = transition_data['FRET_before'] < transition_data['FRET_after'] 
    consecutive_identified = determine_first_transition_in_sequence(transition_data)
    col.append(consecutive_identified)
consecutive_data = pd.concat(col)

consecutive_trans, nonconsecutive_trans, percent_trans_meet_criteria_df = concat_trans_proportion(
    consecutive_data, 
    calculated_transitions_df, 
    FRET_before, 
    FRET_after)

####### calling functions to plot synchronised transition
consecutive_from_dnak_release = filt_df_to_plot(consecutive_trans, FRET_before, FRET_after,'low_to_high', frames_to_plot)
nonconsecutive_from_dnak_release = filt_df_to_plot(nonconsecutive_trans, FRET_before, FRET_after,'low_to_high', frames_to_plot)

plot_synchronised_transition(calculated_transitions_df, consecutive_from_dnak_release, exposure, list_to_drop, order, frames_to_plot, 'consecutive_transitions')
plot_synchronised_transition(calculated_transitions_df, nonconsecutive_from_dnak_release, exposure, list_to_drop, order, frames_to_plot, 'non-consecutive_transition')


########### Code to plot the percentage transition data 

melted_data = percent_trans_meet_criteria_df.melt(id_vars = 'treatment')
fig, ax = plt.subplots()
sns.set_style('ticks',{'grid.linestyle':'--', 'font_scale': 1.5})
sns.barplot(data = melted_data, y = 'value', x = 'variable', hue = 'treatment', palette = 'BuPu', hue_order = order)
plt.xlabel('')
plt.xticks(rotation = 45)
plt.legend(title = '')
plt.ylabel('Proportion of transitions (%)')
fig.savefig(f'{plot_export}/consecutive_transition_summary.svg', dpi = 600)
plt.show()


###############
###############
###############


def ratio_consecutive_to_nonconsecutive(calculated_transitions_df, frames_to_plot):
    consecutive_from_dnak_release = filt_df_to_plot(consecutive_trans, FRET_before, FRET_after,'low_to_high', frames_to_plot)
    nonconsecutive_from_dnak_release = filt_df_to_plot(nonconsecutive_trans, FRET_before, FRET_after,'low_to_high', frames_to_plot)

    test = calculated_transitions_df.iloc[consecutive_from_dnak_release].groupby('treatment_name')['molecule_number'].nunique()/(calculated_transitions_df.iloc[nonconsecutive_from_dnak_release].groupby('treatment_name')['molecule_number'].nunique())
    testies = pd.DataFrame(test).reset_index()
    testies.columns = ['treatment', 'prop_consecutive_dnaK_release']
    return testies


def prop_DnaK_release_events_are_consecutive(calculated_transitions_df, frames_to_plot):
    consecutive_from_dnak_release = filt_df_to_plot(consecutive_trans, FRET_before, FRET_after,'low_to_high', frames_to_plot)
    nonconsecutive_from_dnak_release = filt_df_to_plot(nonconsecutive_trans, FRET_before, FRET_after,'low_to_high', frames_to_plot)
    test = (calculated_transitions_df.iloc[consecutive_from_dnak_release].groupby('treatment_name')['molecule_number'].nunique())/((calculated_transitions_df.iloc[nonconsecutive_from_dnak_release].groupby('treatment_name')['molecule_number'].nunique())+calculated_transitions_df.iloc[consecutive_from_dnak_release].groupby('treatment_name')['molecule_number'].nunique())
    testies = pd.DataFrame(test).reset_index()
    testies.columns = ['treatment', 'prop_consecutive_dnaK_release']
    return testies


def plot_consec_DnaK_release_with_filter(dataframe, datatype = 'Proportion'):
    helpplease = []
    for x, df in enumerate(range(0, 401)):
        if datatype == 'Proportion':
            dfs = prop_DnaK_release_events_are_consecutive(dataframe, x)
        else:
            dfs = ratio_consecutive_to_nonconsecutive(dataframe, x)
        dfs['frames_to_thresh'] = x
        helpplease.append(dfs)
        helpplease_df = pd.concat(helpplease)
        # melty = helpplease_df.melt(id_vars = 'treatment')
    sns.lineplot(data = helpplease_df, x = 'frames_to_thresh', y = 'prop_consecutive_dnaK_release', hue = 'treatment', palette = 'BuPu')
    plt.xlabel('Threshold prior to DnaK release (frames)')
    plt.ylabel(f'{datatype} of transitions (consecutive:non-consecutive)')
    plt.legend(title = '')
    plt.savefig(f'{plot_export}/consecutive_transition_over_frame_threshold_{datatype}.svg', dpi = 600)
    plt.show()
    return

plot_consec_DnaK_release_with_filter(calculated_transitions_df)
