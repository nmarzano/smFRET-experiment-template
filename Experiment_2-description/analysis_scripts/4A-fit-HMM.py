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



output_folder = "Experiment_2-description/python_results/"  ### Change for each experiment
saved_folders = f'{output_folder}/organized_csvs/'
if not os.path.exists(saved_folders):
    os.makedirs(saved_folders)


dict_to_concat = {
    'RNA09-AF488':'Experiment_2-description/raw_data/Chan1_RNA09-AF488_20nM_steady_BC_all',
    # 'RNA09-duplex':'Experiment_1-RNA09-AF488/raw_data/RNA09-DNA02-duplex_pprfret',
    # 'RNA01-AF488':'Experiment_1-RNA09-AF488/raw_data/RNA01-AF488_20nM_pprfret',
    'RNA01-AF488_col':'Experiment_2-description/raw_data/RNA01-AF488_combined-with-230214_pprfret',

 }

def file_reader_3colour(input_folder, data_type):
    """will import data

    Args:
        input_folder (directory): where data is stored
        data_type (str): what data will be used to plot, needs to be either 'hist', 'TDP', 'transition_frequency'
        or 'other'. 

    Returns:
        dataframe: dataframe with data to be used in subseqeunt codes
    """
    if data_type == 'hist':
        filenames = glob.glob(input_folder + "/*.txt")
        dfs = []
        for i, filename in enumerate(filenames):
            molecule_number = re.split(r'(\d+)', filename)[-4]
            movie_number = re.split(r'(\d+)', filename)[-8]
            cum_mol = i+1
            hist_data = pd.read_table(filename, sep="\s+")
            hist_data['molecule number'] = molecule_number
            hist_data['movie number'] = movie_number
            hist_data['cumulative_molecule'] = cum_mol
            dfs.append(hist_data)
        test = pd.concat(dfs)
        test.dropna(axis=1, how='all', inplace = True)
        test.columns = ['Time at 532', 'Frame at 532', 'AF488 at 532', 'Cy3 at 532', 'AF647 at 532','Time at 488', 'Frame at 488', 'AF488 at 488', 'Cy3 at 488', 'AF647 at 488','Time at 488_2', 'Frame at 488_2','FRET AF488 to Cy3', 'Idealized FRET AF488 to Cy3', 'FRET AF488 to AF647', 'Idealized FRET AF488 to AF647', 'Time at 532_2', 'Frame at 532_2', 'FRET Cy3 to AF647', 'Idealized FRET Cy3 to AF647','molecule number', 'movie_number', 'cumulative_molecule']
        test.drop(['Time at 532_2', 'Frame at 532_2', 'Time at 488_2', 'Frame at 488_2'], axis = 1, inplace = True)
        # test['FRET at 488'] = test['Acceptor at 488']/(test['Acceptor at 488']+test['Donor at 488'])
        # test['IntensitySum'] = 
        return pd.DataFrame(test)
    else:
        print('invalid data_type, please set data_type as "hist", "TDP","transition_frequency" or "other" if using for violin or heatmap plots')

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


FRET_bound = compiled_df[compiled_df['FRET Cy3 to AF647']>0.5]
compiled_df['bound'] = np.where(compiled_df['FRET Cy3 to AF647']>0.5, 'RNA-bound', 'Unbound')
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


# for (treatment, molecule), df in compiled_df_HMM.groupby(['treatment','cumulative_molecule']):
#     fig, ax = plt.subplots(3, 1, sharex=True)
#     sns.set_style("whitegrid",{'figure.figsize':(8,5.5), 'grid.linestyle':'--', 'font_scale':1.5} )
#     sns.lineplot(data = df, x = 'Time at 532', y = 'Cy3 at 532', ax = ax[0], color = 'green')
#     sns.lineplot(data = df, x = 'Time at 532', y = 'AF647 at 532', ax = ax[0], color = 'purple')
#     sns.lineplot(data = df, x = 'Time at 488', y = 'AF488 at 488', ax = ax[1], color = 'royalblue')
#     sns.lineplot(data = df, x = 'Time at 488', y = 'Cy3 at 488', ax = ax[1], color = 'olivedrab')
#     sns.lineplot(data = df, x = 'Time at 488', y = 'AF647 at 488', ax = ax[1], color = 'orange')
#     # sns.lineplot(data = df, x = 'Time at 488', y = 'Cy3 FRET cascade', ax = ax[2], color = 'grey')
#     # sns.lineplot(data = df, x = 'Time at 488', y = 'AF647 FRET cascade', ax = ax[2], color = 'orange')
#     sns.lineplot(data = df, x = 'Time at 488', y = 'probe_summed_fluorescence', ax = ax[1], color = 'darkgrey')
#     sns.lineplot(data = df, x = 'Time at 532', y = 'FRET Cy3 to AF647', ax = ax[2], color = 'black')
#     sns.lineplot(data = df, x = 'Time at 532', y = 'e_pred_global', ax = ax[2], color = 'orange')

#     ax[2].set_ylim(0, 1)
#     ax[2].set_xlim(0, 300)
#     for x in ax:
#         [y.set_linewidth(2) for y in x.spines.values()]
#         [y.set_color('black') for y in x.spines.values()]
#     ax[0].set_ylabel('')
#     ax[1].set_ylabel('')
#     ax[2].set_ylabel('FRET')
#     ax[2].set_xlabel('Time (s)')
#     fig.text(0.04, 0.65, 'Fluorescence intensity (a.u.)', ha='center', va='center', rotation='vertical')
#     ax[0].set_title(f'{treatment} molecule {molecule}')
#     plot_dir = f'{plot_export}/{treatment}/'
#     if not os.path.exists(plot_dir):
#         os.makedirs(plot_dir)
#     fig.savefig(f'{plot_dir}/molecule_{molecule}_from{treatment}.svg', dpi = 600)
# plt.show()

