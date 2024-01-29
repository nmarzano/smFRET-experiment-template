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
import os as os
import math 

from Utilities.Data_analysis import file_reader

input_folder = 'Experiment_1-description/python_results'
output_folder = f'{input_folder}/Heatmaps-and-first-transitions'
if not os.path.exists(output_folder):
    os.makedirs(output_folder)
    
transition_type = 'low_to_high'   ##### select either 'low_to_high' or 'high_to_low' for plotting first specified transition
to_filt = ['TREATMENT'] ######### select what datasets you want to keep in your average heatmap


####### import data
compiled_df = pd.read_csv(f'{output_folder}/compiled_df.csv')
col = pd.read_csv(f'{output_folder}/col.csv')
normalised_data = pd.read_csv(f'{output_folder}/normalised_data.csv')
no_filt = compiled_df['treatment_name'].unique().tolist()


#########
######### plot heatmaps and average heatmaps
#########

def plot_heatmap(df, gridsize, bins_hex, plot_type = 'hex'):
    for treatment, dfs in df.groupby('treatment_name'):
        if plot_type == 'hex':
            plt.rcParams['svg.fonttype'] = 'none'
            sns.set(style="whitegrid")
            g = sns.JointGrid(data = dfs, x='time', y='FRET', xlim = (0,100), ylim = (0, 1))
            g.plot_joint(plt.hexbin, gridsize=(gridsize, gridsize), cmap='ocean_r', mincnt=0, bins=bins_hex)
            g.plot_marginals(sns.histplot, kde=True, bins=20)
            plt.savefig(f'{output_folder}/Heatmap_{treatment}_{plot_type}.svg', dpi = 600)
        if plot_type == 'kde':
            plt.rcParams['svg.fonttype'] = 'none'
            sns.set_style("whitegrid", {'grid.linestyle':'--', 'axes.linewidth':20, 'axes.color':'black', 'axes.edgecolor': 'black', 'font.size':10})
            fig = sns.jointplot(data = dfs, x='time', y='FRET', xlim = (0,200), ylim = (0, 1), alpha = 0.05, color = '#2AA6CF', marginal_kws = dict(bins = 20, kde = True))   
            fig.plot_joint(sns.kdeplot, cmap = 'mako') 
            fig.ax_joint.spines['top'].set_visible(True)
            fig.ax_joint.spines['right'].set_visible(True)
            fig.ax_marg_x.spines['left'].set_visible(True)
            fig.ax_marg_x.spines['top'].set_visible(True)
            fig.ax_marg_x.spines['right'].set_visible(True)
            fig.ax_marg_y.spines['top'].set_visible(True)
            fig.ax_marg_y.spines['right'].set_visible(True)
            fig.ax_marg_y.spines['bottom'].set_visible(True)
            plt.savefig(f'{output_folder}/Heatmap_{treatment}_{plot_type}.svg', dpi = 600)
        plt.show()
    return


def plot_average_FRET_over_time(df, filt, ci = 'sd', x_axis = 'time', subplot = False):
    sns.set(style = "ticks", font_scale = 1)
    filt_data = df[df['treatment_name'].isin(filt)]
    if x_axis == 'time':
        fig, ax = plt.subplots()
        sns.lineplot(data=filt_data, x = x_axis, y = 'FRET', ci = ci, hue = 'treatment_name')
        plt.xlim(0, 200, 10)
        ax.axvline(x = 10, linestyle = '-', color = 'grey')
        ax.axvspan(0, 10, facecolor='grey', alpha = .2)   
        plt.ylim(0, 1, 10)
        plt.xlabel('Time (s)')
        plt.legend(title = '', fontsize = 'small')
        plt.savefig(f'{output_folder}/Average_Heatmap_{filt}.svg', dpi = 600)
    if (x_axis == 'normalised_to_event' and subplot == True):
        nrow = math.ceil(len(filt)/2)
        fig, axes = plt.subplots(nrow, 2, figsize = (8, 3*nrow), sharex = True, sharey = True)
        axes = axes.flatten()
        for i, treatment in enumerate(filt):
            sns.lineplot(data=df[df['treatment_name'] == treatment], x = x_axis, y = 'FRET', ci = ci, hue = 'treatment_name', ax = axes[i])
            axes[i].axvline(x = 0, linestyle = '--', color = 'grey')
            plt.xlim(-30, 30, 10)
            plt.ylim(0, 1, 10)
            axes[i].set_xlabel('Time (s)')
            axes[i].set_title(f"{treatment} (n = {str(normalised_data[normalised_data['treatment_name'] == treatment]['molecule_number'].nunique())})")
            axes[i].get_legend().remove()
        if len(filt) < len(axes):
            axes[-1].remove()
            axes[3].set_xlabel('Time (s)')
        plt.savefig(f'{output_folder}/Traces_normalised_to_first_trans.svg', dpi = 600)
        plt.show()
    elif x_axis == 'normalised_to_event':
        for treatment, df in filt_data.groupby('treatment_name'):
            sns.lineplot(data=df[df['treatment_name'] == treatment], x = x_axis, y = 'FRET', ci = ci, hue = 'treatment_name')
            plt.xlim(-30, 30, 10)
            plt.axvline(x = 0, linestyle = '--', color = 'grey')
            plt.ylim(0, 1, 10)
            plt.xlabel('Time (s)')
            plt.legend(title = '', fontsize = 'small', loc = 'upper left')
            plt.savefig(f'{output_folder}/Traces_normalised_to_first_trans_{treatment}.svg', dpi = 600)
            plt.show()
    plt.show()
    return



def plot_first_specified_transition(df, trans_type):
    plot1 = plt.figure()
    sns.set(style = 'ticks', font_scale = 1)
    sns.violinplot(data = df, y = 'cum_sum', x = 'treatment', scale = 'width')
    sns.stripplot(data = df, y = 'cum_sum', x = 'treatment', color ='black', alpha = 0.5)
    plt.ylabel('Time for RNA binding (s)')
    plt.xlabel('')
    plot1.savefig(f'{output_folder}/time_until_first_{trans_type}_transition.svg', dpi = 600)
    plt.show()


plot_heatmap(compiled_df, 100, 80, 'kde')
plot_first_specified_transition(col, transition_type)
plot_average_FRET_over_time(normalised_data, no_filt, 'sd', 'time') ##### change 'no_filt' to 'to_filt' to include only datasets that were mentioned above
plot_average_FRET_over_time(normalised_data, no_filt, 'sd', 'normalised_to_event', False)  #### replace with True to plot subplots

print(col.groupby('treatment')['cum_sum'].mean())
print(col.groupby('treatment')['cum_sum'].sem())






