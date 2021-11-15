from scipy.stats import norm
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import pandas as pd
import numpy as np
import glob
import os

output_folder = 'Experiment_1-description/python_results'
data_paths = {
    "key":"Experiment_1-description/raw_data/_PATH_D051521_T1159.dat",
    'key':'Experiment_1-description/raw_data/PATH_D051521_T1447.dat',
}


def remove_outliers(compiled_TDP):
    outliers = compiled_TDP[(compiled_TDP["FRET_before"] < -0.5)|(compiled_TDP["FRET_before"] > 1.5)|(compiled_TDP["FRET_after"] < -0.5) | (compiled_TDP["FRET_after"] > 1.5)].index
    compiled_TDP.drop(outliers, inplace = True)
    return compiled_TDP

def calculate_dwell_time(df):
    df_test2 = []
    for Molecule, dfs in df.groupby('Molecule'):
        frame_length = len(dfs)
        test = dfs.groupby([dfs['Idealized_FRET'].ne(dfs['Idealized_FRET'].shift()).cumsum(), 'Idealized_FRET']).size()
        test = test.reset_index(level = 1, drop = False)
        test['Molecule'] = Molecule
        test['number_of_frames'] = frame_length
        df_test2.append(test)
    df_test3 = pd.concat(df_test2)
    df_test3.columns = ['FRET_state', 'Time', 'Molecule', 'number_of_frames']
    df_test3 = df_test3.reset_index().drop('Idealized_FRET', axis = 1)
    df_test4 = df_test3[df_test3.groupby('Molecule').Molecule.transform('count') > 1]
    return df_test4

def generate_transitions(df):
    df_toconcat = []
    for molecule, dfs in df.groupby('Molecule'):
        poop = dfs.assign(FRET_after = dfs.FRET_state.shift(-1)).dropna()
        df_toconcat.append(poop)
    compiled_df = pd.concat(df_toconcat).reset_index(drop = True)
    compiled_final = compiled_df[['Molecule', 'FRET_state', 'FRET_after', 'Time', 'number_of_frames']]
    compiled_final.columns = ['Molecule', 'FRET_before', 'FRET_after', 'Time', 'number_of_frames']
    return compiled_final

from Utilities.Data_analysis import filter_TDP, file_reader, count_filt_mol

compiled_TDP = []
for data_name, data_path in data_paths.items():
    imported_TDPdata = file_reader(data_path, 'TDP')
    df1 = calculate_dwell_time(imported_TDPdata)
    df2 = generate_transitions(df1)
    cleaned_TDPdata = remove_outliers(df2)
    cleaned_TDPdata["treatment_name"] = data_name
    compiled_TDP.append(cleaned_TDPdata)
compiled_TDP = pd.concat(compiled_TDP)   
compiled_TDP.to_csv(f'{output_folder}/TDP_cleaned.csv', index = False)

raw_TDPdata = []
for data_name, data_path in data_paths.items():
    imported_TDPdata = file_reader(data_path, 'TDP')
    df1 = calculate_dwell_time(imported_TDPdata)
    df2 = generate_transitions(df1)
    # cleaned_TDPdata = remove_outliers(df2)
    df2["treatment_name"] = data_name
    raw_TDPdata.append(df2)
compiled_raw_TDP = pd.concat(raw_TDPdata)   
compiled_raw_TDP.to_csv(f'{output_folder}/TDP_raw.csv', index = False)

filtered_data = filter_TDP(compiled_TDP, 0.3)  ##### number refers to the FRET threshold to filter data. Will only include molecules that go below this set threshold. Will default to 0.3 
filtered_data.to_csv(f'{output_folder}/TDP_cleaned_filt.csv', index = False)

############### calculate mol count
FRET_value = 0.6
molcount = count_filt_mol(compiled_TDP, FRET_value, data_paths, 1)  ## last number refers to the index of the treatment you want to subtract from others
molcount.to_csv(f'{output_folder}/mol_below_{FRET_value}.csv', index = None)


#### .rename(columns = {1:"test", 3:"test2"}) ## can rename individually if needed
# compiled_TDP.columns = ['Molecule', 'FRET before transition', 'FRET after transition', 'Time', "treatment_name"]

