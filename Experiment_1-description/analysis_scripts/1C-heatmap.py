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

output_folder = 'Experiment_1-description/python_results'
exposure = 0.200  # in seconds
frame_rate = 1/exposure
FRET_thresh = 0.6 #### Used to filt the threshold for transitions


data_paths = {
    'treatment':('treatment_name', 'Experiment_1-description/raw_data/treatment_to_plot'),
}

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
            return compiled[(compiled['FRET'] > -0.5) & (compiled['FRET'] < 1.5)].copy()
        if data_type == "idealized":
            return compiled[(compiled[4] > -0.5) & (compiled[4] < 1.5)].copy()
    elif plot_type == 'TDP':
        outliers = compiled[(compiled["FRET before transition"] < -0.5)|(compiled["FRET before transition"] > 1.5)|(compiled["FRET after transition"] < -0.5) | (compiled["FRET after transition"] > 1.5)].index
        compiled.drop(outliers, inplace = True)
        return compiled
    else:
        print('invalid plot type, please set plot_type as "hist" or "TDP" - you idiot')

def plot_heatmap(df, gridsize, bins_hex):
    for treatment, dfs in df.groupby('treatment_name'):
        plt.rcParams['svg.fonttype'] = 'none'
        sns.set(style="whitegrid")
        g = sns.JointGrid(data = dfs, x='time', y='FRET', xlim = (0,100), ylim = (0, 1))
        g.plot_joint(plt.hexbin, gridsize=(gridsize, gridsize), cmap='ocean_r', mincnt=0, bins=bins_hex)
        g.plot_marginals(sns.histplot, kde=True, bins=20)
        identifier = treatment
        plt.savefig(f'{output_folder}/Heatmap_{identifier}.svg', dpi = 600)
        plt.show()
    return

compiled_df = []
for data_name, (label, data_path) in data_paths.items():
    imported_data = file_reader(data_path, 'heatmap', frame_rate)
    cleaned_raw = remove_outliers(imported_data, 'hist')    #### add "idealized" after imported_data to get idealized histograms
    cleaned_raw["treatment_name"] = data_name
    compiled_df.append(cleaned_raw)
compiled_df = pd.concat(compiled_df)   #### .rename(columns = {1:"test", 3:"test2"}) ## can rename individually if needed
compiled_df.columns = ["frames", "donor", "acceptor", "FRET", "idealized FRET", 'molecule_number', 'time', "treatment_name"]
plot_heatmap(compiled_df, 100, 80)




def calculate_dwell_time(df):
    """Function to convert raw idealized data to a form in which the duration of each idealized state is calculated

    Args:
        df (dataframe): dataframe containing each molecule and the idealized fret 

    Returns:
        dataframe: returns dataframe containing the duration of each FRET state from each molecule
    """
    df_test2 = []
    for Molecule, dfs in df.groupby('molecule_number'):
        frame_length = len(dfs)
        test = dfs.groupby([dfs['idealized FRET'].ne(dfs['idealized FRET'].shift()).cumsum(), 'idealized FRET']).size()
        test = test.reset_index(level = 1, drop = False)
        test['Molecule'] = Molecule
        test['number_of_frames'] = frame_length
        df_test2.append(test)
    df_test3 = pd.concat(df_test2)
    df_test3.columns = ['FRET_state', 'Time', 'Molecule', 'number_of_frames']
    df_test3 = df_test3.reset_index().drop('idealized FRET', axis = 1)
    return df_test3[df_test3.groupby('Molecule').Molecule.transform('count') > 1]


def generate_transitions(df):
    """Converts the duration of each FRET state into a transition, whereby the FRET state before, the FRET state after
    and the duration of the FRET state before a transition is given in a single line. Each line represents a single transition.

    Args:
        df (dataframe): dataframe generated following 'calculate_dwell_time' function in which the duration of a certain
        FRET state is given for each FRET state for all molecules

    Returns:
        dataframe: returns a dataframe in which each row represents a transition, with FRET before transition, FRET after transition
        and duration of FRET state before transition (given in number of frames in column Time) provided
    """
    df_toconcat = []
    for molecule, dfs in df.groupby('Molecule'):
        thing1 = dfs.assign(FRET_after = dfs.FRET_state.shift(-1)).dropna()
        df_toconcat.append(thing1)
    compiled_df = pd.concat(df_toconcat).reset_index(drop = True)
    compiled_final = compiled_df[['Molecule', 'FRET_state', 'FRET_after', 'Time', 'number_of_frames']]
    compiled_final.columns = ['Molecule', 'FRET_before', 'FRET_after', 'Time', 'number_of_frames']
    return compiled_final

def remove_outliers(compiled_TDP):
    outliers = compiled_TDP[(compiled_TDP["FRET_before"] < -0.5)|(compiled_TDP["FRET_before"] > 1.5)|(compiled_TDP["FRET_after"] < -0.5) | (compiled_TDP["FRET_after"] > 1.5)].index
    compiled_TDP.drop(outliers, inplace = True)
    return compiled_TDP


def filter_FRET_trans_if(dfs, thresh, trans_type = 'low_to_high'):
    comb = []
    for Molecule, df in dfs.groupby('Molecule'):
        df['cum_sum'] = df['time (s)'].cumsum()
        comb.append(df)
    combined = pd.concat(comb)
    if trans_type == 'high_to_low':
        filt_data = combined[(combined['FRET_before'] > thresh) & (combined['FRET_after'] < thresh)]
    elif trans_type == 'low_to_high':
        filt_data = combined[(combined['FRET_before'] < thresh) & (combined['FRET_after'] > thresh)]
    return filt_data




def select_first_transition(dfs):
    first_trans = []
    for molecule, df in dfs.groupby('Molecule'):
        first_trans_above_thresh = df[df['cum_sum'] == df['cum_sum'].min()]
        first_trans.append(first_trans_above_thresh)
    return pd.concat(first_trans)


compiled_filt = []
for treatment, df in compiled_df.groupby('treatment_name'):
    treatment_df = compiled_df[compiled_df['treatment_name'] == treatment]
    treatment_df2 = treatment_df.filter(items = ['idealized FRET','molecule_number'])
    treatment_df3 = calculate_dwell_time(treatment_df2)
    treatment_transitions = generate_transitions(treatment_df3)
    treatment_cleaned_transitions = remove_outliers(treatment_transitions)
    treatment_cleaned_transitions['time (s)'] = treatment_cleaned_transitions['Time'] * exposure
    treatment_cumsum = filter_FRET_trans_if(treatment_cleaned_transitions, FRET_thresh) ##### add 'high_to_low' to look at how long it takes for high-low transitions occur
    treatment_first_transition = select_first_transition(treatment_cumsum)
    treatment_first_transition['treatment'] = treatment
    compiled_filt.append(treatment_first_transition)
col = pd.concat(compiled_filt)

plot1 = plt.figure()
sns.set_style('ticks')
sns.violinplot(data = col, y = 'cum_sum', x = 'treatment', scale = 'width')
sns.stripplot(data = col, y = 'cum_sum', x = 'treatment', color ='black', alpha = 0.5)
plt.ylabel('Time for RNA binding (s)')
plot1.savefig(f'{output_folder}/time_until_first_transition_above_thresh.svg', dpi = 600)


col.groupby('treatment')['cum_sum'].mean()
col.groupby('treatment')['cum_sum'].sem()

