from scipy.stats import norm
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import pandas as pd
import numpy as np
import shutil
import os
import glob
import re
from deepfret_nm import DeepFRET_NM
from Utilities.Data_analysis import file_reader_3colour

output_folder = "Experiment_1-description/python_results/"  ### Change for each experiment
saved_folders = f'{output_folder}/organized_csvs/'
if not os.path.exists(saved_folders):
    os.makedirs(saved_folders)

FRET_thresh = 0.5

dict_to_concat = {
    'treatment1':'Experiment_1-description/raw_data/treatment1_data',
    'treatment2':'Experiment_1-description/raw_data/treatment2_data',

 }

compiled_df = []
for data_name, data_path in dict_to_concat.items():
    imported_data = file_reader_3colour(data_path, 'hist')
    imported_data["treatment"] = data_name
    compiled_df.append(imported_data)
compiled_df = pd.concat(compiled_df) 
compiled_df['Cy3 FRET'] = compiled_df['Cy3 at 488']/(compiled_df['Cy3 at 488']+compiled_df['AF488 at 488']+compiled_df['AF647 at 488'])
compiled_df['AF647 FRET'] = compiled_df['AF647 at 488']/(compiled_df['Cy3 at 488']+compiled_df['AF488 at 488']+compiled_df['AF647 at 488'])
compiled_df['Cy3 FRET cascade'] = (compiled_df['AF647 at 488'] + compiled_df['Cy3 at 488'])/(compiled_df['Cy3 at 488']+compiled_df['AF488 at 488']+compiled_df['AF647 at 488'])
compiled_df['AF647 FRET cascade'] = (compiled_df['AF647 at 488'])/(compiled_df['Cy3 at 488']+compiled_df['AF647 at 488'])
compiled_df['Cy3 FRET Lee'] = (compiled_df['Cy3 at 488'])/(compiled_df['Cy3 at 488']+(compiled_df['AF488 at 488']*(1 - compiled_df['FRET Cy3 to AF647'])))
compiled_df['probe_summed_fluorescence'] = compiled_df['Cy3 at 488'] + compiled_df['AF488 at 488'] + compiled_df['AF647 at 488']


FRET_bound = compiled_df[compiled_df['FRET Cy3 to AF647']>FRET_thresh]
compiled_df['bound'] = np.where(compiled_df['FRET Cy3 to AF647']>FRET_thresh, 'RNA-bound', 'Unbound')
print(FRET_bound.groupby('treatment')['cumulative_molecule'].nunique())
print(compiled_df.groupby('treatment')['cumulative_molecule'].nunique())
compiled_df_filt = compiled_df.drop(['Frame at 532'], axis = 1)



list(compiled_df_filt.columns)
compiled_df_filt = compiled_df_filt[[
 'Time at 532',
 'Cy3 at 532',
 'AF647 at 532',
 'AF488 at 532',
 'Time at 488',
 'Frame at 488',
 'AF488 at 488',
 'Cy3 at 488',
 'AF647 at 488',
 'FRET AF488 to Cy3',
 'Idealized FRET AF488 to Cy3',
 'FRET AF488 to AF647',
 'Idealized FRET AF488 to AF647',
 'FRET Cy3 to AF647',
 'Idealized FRET Cy3 to AF647',
 'molecule number',
 'movie_number',
 'cumulative_molecule',
 'treatment',
 'Cy3 FRET',
 'AF647 FRET',
 'Cy3 FRET cascade',
 'AF647 FRET cascade',
 'Cy3 FRET Lee',
 'probe_summed_fluorescence',
 'bound'
 ]]


for (treatment, molecule), df in compiled_df_filt.groupby(['treatment', 'cumulative_molecule']):
    df_filt = df[['Time at 532', 'Cy3 at 532', 'AF647 at 532']]
    df_filt.to_csv(f'{saved_folders}/{treatment}_molecule{molecule}.dat', index = False,header = False, sep='t')


compiled_HMM = []
for (treatment, molecule), df in compiled_df_filt.groupby(['treatment', 'cumulative_molecule']):
    df
    trace = DeepFRET_NM.import_data(f'Experiment_1-RNA09-AF488/raw_data/organized_csvs/{treatment}_molecule{molecule}.dat')
    traces = [trace]
    df['e_pred_global'] = DeepFRET_NM.fit_HMM_NM(traces)['e_pred_global']
    compiled_HMM.append(df)
compiled_df_HMM = pd.concat(compiled_HMM)
compiled_df_HMM.to_csv(f'{output_folder}/compiled_df_HMM.csv', index=False)
print(compiled_df_HMM)

