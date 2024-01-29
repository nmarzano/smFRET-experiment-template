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
plot_export = f'{output_folder}/histograms/'
if not os.path.exists(plot_export):
    os.makedirs(plot_export)

FRET_thresh = 0.5

###### import data
compiled_df_HMM = pd.read_csv(f'{output_folder}/compiled_df_HMM.csv')

###### Identify when molecule is in a certain state (e.g., bound or unbound) based on FRET efficiency and label when it is in that state.
FRET_bound = compiled_df_HMM[compiled_df_HMM['FRET Cy3 to AF647']>FRET_thresh]
compiled_df_HMM['bound'] = np.where(compiled_df_HMM['FRET Cy3 to AF647']>FRET_thresh, 'RNA-bound', 'Unbound')

colors = {
    "treatment1": "black", 
    "treatment2":"grey", 
}
###########
########### plots two graphs, one with the PPR FRET and another with whatever FRET you want to look at
###########

def plot_FRET_multiple(df, top_FRET, bottom_FRET):
    data_hist = df[(df[top_FRET] > -0.2) & (df[top_FRET] < 1.2) & (df[bottom_FRET] > -0.2) & (df[bottom_FRET] < 1.2)]
    fig, axes = plt.subplots(2, 1, sharex =  True)
    sns.set_style("ticks",{'figure.figsize':(8,5.5),'font_scale':1.5} )
    sns.histplot(data = data_hist, x = top_FRET, binrange = (0, 1), binwidth = 0.05, kde = True, ax = axes[0], stat = 'density', fill = False, hue = 'treatment', common_norm = False, palette = 'gray')
    sns.histplot(data = data_hist, x = bottom_FRET, binrange = (0, 1), binwidth = 0.05, kde = True, ax = axes[1], stat = 'density', fill = False, hue = 'treatment', common_norm = False, palette = 'gray')
    for x in axes:
        [y.set_linewidth(1) for y in x.spines.values()]
        [y.set_color('black') for y in x.spines.values()]
    axes[0].set_title(top_FRET)
    axes[1].set_title(bottom_FRET)
    plt.xlabel('FRET')
    plt.xlim(0, 1)
    plt.savefig(f'{plot_export}/FRET-histograms.svg', dpi = 600)
    plt.show()


plot_FRET_multiple(compiled_df_HMM, 'FRET Cy3 to AF647', 'FRET Cy3 to AF647')
plot_FRET_multiple(FRET_bound, 'Cy3 FRET cascade', 'AF647 FRET cascade')

################
################ Plots FRET distribution. 'Bar' is to plot multiple treatments in the one graph
################ 'anything else' is to plot each treatment independently. 
################

def plot_hist_type(df, kind = 'kde', palette = colors):
    plot_hist, ax = plt.subplots()
    sns.set_style("ticks",{'font_scale':1} )
    plt.xlim(0, 1, 10)
    plt.xlabel("FRET")
    [x.set_linewidth(2) for x in ax.spines.values()]
    [x.set_color('black') for x in ax.spines.values()]
    if kind == 'kde':
        sns.kdeplot(
            data = df, 
            palette = palette, 
            x = "FRET Cy3 to AF647",
            hue="treatment",
            common_norm=False, 
            fill = True, 
            linewidth = 1.5, 
            alpha = 0.25)
        plt.show()
    if kind == 'bar':
        sns.histplot(
            data = df, 
            x = "FRET Cy3 to AF647",
            common_norm=False, 
            stat = 'density',
            hue = 'treatment',
            palette = palette,
            binwidth = 0.05,
            fill = False, 
            kde = True,
            linewidth = 1.5, 
            alpha = 0.25)
        plt.show()
    else: 
        for treatment, dfs in df.groupby('treatment'):
            df_filt = dfs[dfs['treatment']==treatment]
            sns.histplot(
                data = df_filt, 
                x = "FRET Cy3 to AF647",
                common_norm=False, 
                stat = 'density',
                color = 'black',
                binwidth = 0.05,
                fill = False, 
                kde = True,
                linewidth = 1.5, 
                alpha = 0.25)
            plt.xlim(0, 1, 10)
            plt.xlabel("FRET")
            [x.set_linewidth(2) for x in ax.spines.values()]
            [x.set_color('black') for x in ax.spines.values()]
            plot_hist.savefig(f'{plot_export}/Histogram_{treatment}.svg', dpi = 600)
            plt.show()
    plot_hist.savefig(f'{plot_export}/Histogram_{kind}.svg', dpi = 600)
    


plot_hist_type(compiled_df_HMM, 'bar')


###############
############### plots the fluorescence intensity of a single dye for each oligo
############### in the PPR bound or unbound state
###############


palette_dict_treatment = {
    'RNA01-AF488':{'RNA-bound':'black', 'Unbound':'grey'}, 
    'RNA09-AF488':{'RNA-bound':'royalblue', 'Unbound':'skyblue'}, 
    'RNA09-duplex':{'RNA-bound':'purple', 'Unbound':'orchid'}, 
    'RNA01-AF488_col':{'RNA-bound':'black', 'Unbound':'grey'}, 


    }

def plot_intensity_for_treatments(df, intensity_type, palette):
    treatment_list = list(df['treatment'].unique())
    fig, axes = plt.subplots(len(treatment_list), 1, sharex=True)
    sns.set_style('ticks',{'grid.linestyle':'--', 'font_scale': 1.5})
    for (i, label) in enumerate(treatment_list):
        treatment_data = df[df['treatment'] == label]
        sns.histplot(data = treatment_data, 
                     x = intensity_type, 
                     kde = True, 
                     stat = 'density', 
                     fill = False, 
                     common_norm = False, 
                     palette = palette[label], 
                     binwidth = 250,
                     binrange = (-5000, 20000),
                     hue = 'bound',
                     ax = axes[i])
        axes[i].set_title(f'{label}')
        axes[i].legend(['Unbound', 'RNA-bound'])   
    axes[1].set_xlabel(intensity_type)
    plt.tight_layout()
    plt.xlim(-5000, 20000)
    plt.show()
    return fig, axes

fig, axes = plot_intensity_for_treatments(compiled_df_HMM, 'AF647 at 488', palette_dict_treatment)



#############
############# For each treatment it will plot the fluorescence intensity of the 
############# dyes excited by 488 nm laser in the bound and unbound PPR state
#############

flurophore_palette = {
    'AF488 at 488':{'RNA-bound':'#0077b6', 'Unbound':'#caf0f8'}, 
    'Cy3 at 488':{'RNA-bound':'#31572c', 'Unbound':'#dde5b6'}, 
    'AF647 at 488':{'RNA-bound':'#ae2012', 'Unbound':'#e9d8a6'}, 
    # 'probe_summed_fluorescence':{'RNA-bound':'#ae2012', 'Unbound':'#e9d8a6'}, 

    }


def plot_each_intensity(df, treatment):
    treatment_list = ['AF488 at 488', 'Cy3 at 488', 'AF647 at 488']
    fig, axes = plt.subplots(3, 1, sharex=True)
    sns.set_style('ticks',{'grid.linestyle':'--', 'font_scale': 1.5})
    for (i, label) in enumerate(treatment_list):
        sns.histplot(data = df, 
                     x = label, 
                     kde = True, 
                     stat = 'density', 
                     fill = False, 
                     common_norm = False, 
                     palette = flurophore_palette[label], 
                     binwidth = 500,
                     binrange = (-5000, 20000),
                     hue = 'bound',
                     hue_order = ['RNA-bound', 'Unbound'],
                     ax = axes[i], 
                     legend = True)
        axes[i].set_title(f'{label} for {treatment}')
    axes[2].set_xlabel('Total fluorescence (a.u.)')
    plt.tight_layout()
    plt.xlim(-5000, 20000)
    fig.savefig(f'{plot_export}/intensity_bound_{treatment}.svg', dpi = 600)
    plt.show()
    return fig, axes


for treatment, df in compiled_df_HMM.groupby('treatment'):
    fig, axes = plot_each_intensity(df, treatment)
