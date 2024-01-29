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



output_folder = "Experiment_1-description/python_results/"  ### Change for each experiment
plot_export = f'{output_folder}/traces/'
if not os.path.exists(plot_export):
    os.makedirs(plot_export)

compiled_df_HMM = pd.read_csv(f'{output_folder}/compiled_df_HMM.csv')

for (treatment, molecule), df in compiled_df_HMM.groupby(['treatment','cumulative_molecule']):
    fig, ax = plt.subplots(3, 1, sharex=True)
    sns.set_style("whitegrid",{'figure.figsize':(8,5.5), 'grid.linestyle':'--', 'font_scale':1.5} )
    sns.lineplot(data = df, x = 'Time at 532', y = 'Cy3 at 532', ax = ax[0], color = 'green')
    sns.lineplot(data = df, x = 'Time at 532', y = 'AF647 at 532', ax = ax[0], color = 'purple')
    sns.lineplot(data = df, x = 'Time at 488', y = 'AF488 at 488', ax = ax[1], color = 'royalblue')
    sns.lineplot(data = df, x = 'Time at 488', y = 'Cy3 at 488', ax = ax[1], color = 'olivedrab')
    sns.lineplot(data = df, x = 'Time at 488', y = 'AF647 at 488', ax = ax[1], color = 'orange')
    # sns.lineplot(data = df, x = 'Time at 488', y = 'Cy3 FRET cascade', ax = ax[2], color = 'grey')
    # sns.lineplot(data = df, x = 'Time at 488', y = 'AF647 FRET cascade', ax = ax[2], color = 'orange')
    sns.lineplot(data = df, x = 'Time at 488', y = 'probe_summed_fluorescence', ax = ax[1], color = 'darkgrey')
    sns.lineplot(data = df, x = 'Time at 532', y = 'FRET Cy3 to AF647', ax = ax[2], color = 'black')
    sns.lineplot(data = df, x = 'Time at 532', y = 'e_pred_global', ax = ax[2], color = 'orange')

    ax[2].set_ylim(0, 1)
    ax[2].set_xlim(0, 300)
    for x in ax:
        [y.set_linewidth(2) for y in x.spines.values()]
        [y.set_color('black') for y in x.spines.values()]
    ax[0].set_ylabel('')
    ax[1].set_ylabel('')
    ax[2].set_ylabel('FRET')
    ax[2].set_xlabel('Time (s)')
    fig.text(0.04, 0.65, 'Fluorescence intensity (a.u.)', ha='center', va='center', rotation='vertical')
    ax[0].set_title(f'{treatment} molecule {molecule}')
    plot_dir = f'{plot_export}/{treatment}/'
    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)
    fig.savefig(f'{plot_dir}/molecule_{molecule}_from{treatment}.svg', dpi = 600)
plt.show()

