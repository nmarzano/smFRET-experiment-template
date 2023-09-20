from scipy.stats import norm
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import pandas as pd
import numpy as np
import Utilities.Data_analysis as uda
import os

output_folder = "Experiment_1-description/python_results"
plot_export = f'{output_folder}/Residence_time_plots/'
if not os.path.exists(plot_export):
    os.makedirs(plot_export)
FRET_thresh = 0.2 #### FRET value at which to filter data above or below. 
cumulative_dwell = pd.read_csv(f'{output_folder}/cumulative_dwell.csv')
order = ['', '', '', '', '']

data_paths_violin = {
    f"{order[0]}":f"{output_folder}/Dwell_times/Filtered_dwelltime_{order[0]}.csv",
    f"{order[1]}":f"{output_folder}/Dwell_times/Filtered_dwelltime_{order[1]}.csv",
    f"{order[2]}":f"{output_folder}/Dwell_times/Filtered_dwelltime_{order[2]}.csv",
    f"{order[3]}":f"{output_folder}/Dwell_times/Filtered_dwelltime_{order[3]}.csv",
    f"{order[4]}":f"{output_folder}/Dwell_times/Filtered_dwelltime_{order[4]}.csv",
    # f"{order[5]}":"Figure3b-overhangs/python_results/Dwell_times/Filtered_dwelltime_RNA09.csv",
}

dict_key = list(data_paths_violin.keys())
#### Directory for above should come from the Dwell_times folder in python_results

def compiled(df, data_name):
    """Will filter transitions dependent on a threshold defined above as FRET_thresh to calculate residenc time for each transition class

    Args:
        df (dataframe): dataset containing the residence times  for each treatment
        data_name (_type_): treatment name  

    Returns:
        dataframe: compiles all transition classes (with residence times) from all treatments together
    """
    violin_data_lowtolow = pd.DataFrame(df[f"< {FRET_thresh} to < {FRET_thresh}"])
    violin_data_lowtolow.columns = ["y_axis"]
    violin_data_lowtolow["transition_type"] = f"< {FRET_thresh} to < {FRET_thresh}"
    violin_data_lowtolow["treatment"] = data_name

    violin_data_lowtohigh = pd.DataFrame(df[f"< {FRET_thresh} to > {FRET_thresh}"])
    violin_data_lowtohigh.columns = ["y_axis"]
    violin_data_lowtohigh["transition_type"] = f"< {FRET_thresh} to > {FRET_thresh}"
    violin_data_lowtohigh["treatment"] = data_name

    violin_data_hightohigh = pd.DataFrame(df[f"> {FRET_thresh} to > {FRET_thresh}"])
    violin_data_hightohigh.columns = ["y_axis"]
    violin_data_hightohigh["transition_type"] = f"> {FRET_thresh} to > {FRET_thresh}"
    violin_data_hightohigh["treatment"] = data_name

    violin_data_hightolow = pd.DataFrame(df[f"> {FRET_thresh} to < {FRET_thresh}"])
    violin_data_hightolow.columns = ["y_axis"]
    violin_data_hightolow["transition_type"] = f"> {FRET_thresh} to < {FRET_thresh}"
    violin_data_hightolow["treatment"] = data_name
    return pd.concat(
        [
            violin_data_lowtolow,
            violin_data_lowtohigh,
            violin_data_hightohigh,
            violin_data_hightolow,
        ]
    )


test = []
for data_name, data_path in data_paths_violin.items():
    data = uda.file_reader(data_path, 'other')
    compiled_data = compiled(data, data_name)
    test.append(compiled_data)
final = pd.concat(test)
final["y_axis_log10"] = np.log10(final['y_axis']) ## if need to plot in log scale

############# Code to plot data
############# NEED TO MODIFY BASED ON DATASET but can COPY FROM 1A-PLOT-HISTOGRAM
colors_violin =  {
    'RNA01': 'black', 
    'RNA09':"#ee9b00", 
    'RNA10':"#ca6702", 
    'RNA22':"#ae2012", 
    'RNA22duplex':"royalblue",
    # 'DNA02':"grey"
}

dict_for_label = {
f"< {FRET_thresh} to < {FRET_thresh}":'$T_{low-low}$', 
f"< {FRET_thresh} to > {FRET_thresh}":'$T_{low-high}$',
f"> {FRET_thresh} to > {FRET_thresh}":'$T_{high-high}$',
f"> {FRET_thresh} to < {FRET_thresh}":'$T_{high-low}$'
}


final['transition_name'] = final['transition_type'].map(dict_for_label)


def plot_violin(data, scale = "y_axis"):
    if scale == "y_axis":
        plt.rcParams['svg.fonttype'] = 'none'
        plot1 = plt.figure()
        sns.set(style = "ticks", font_scale = 1)
        sns.violinplot(
            data=data, 
            x= "transition_name", 
            y = "y_axis",
            palette= colors_violin, 
            hue = "treatment", 
            log_scale = True,
            cut = 0)
        plt.ylabel("Residence time (s)")
        plt.xlabel("Transition class")
        plt.legend(title = '',loc = "upper right")
        plot1.savefig(f"{plot_export}/Violin_plot_normal.svg", dpi = 600)
        plt.show()
    if scale == "y_axis_log10":
        plt.rcParams['svg.fonttype'] = 'none'
        plot2 = plt.figure()
        sns.set(style = "ticks", font_scale = 1)
        sns.violinplot(
            data=data, 
            x= "transition_name", 
            y = "y_axis_log10",
            palette= 'mako', 
            hue = "treatment", 
            log_scale = True)
        plt.ylabel("Log residence time (s)")
        plt.xlabel("Transition class")
        plt.legend(title = '',loc = "upper left", bbox_to_anchor = (1,1), ncol =1)
        plot2.savefig(f"{plot_export}/Violin_plot_log.svg", dpi = 600)
        plt.show()
    if scale == 'split':     
        f, (ax_top, ax_bottom) = plt.subplots(ncols=1, nrows=2, sharex=True, gridspec_kw={'hspace':0.05})
        sns.set(style = 'ticks')
        sns.violinplot(x="transition_name", y="y_axis", hue="treatment",data=data, ax=ax_top, palette = 'BuPu', scale = 'width')
        sns.violinplot(x="transition_name", y="y_axis", hue="treatment",data=data, ax=ax_bottom, cut = 0, palette = 'BuPu', scale = 'width')
        ax_top.set_ylim(bottom=40)   # those limits are fake
        ax_bottom.set_ylim(0,40)
        # sns.despine(ax=ax_bottom)
        # sns.despine(ax=ax_top, bottom=True)
        ax = ax_top
        d = .015  # how big to make the diagonal lines in axes coordinates
        # arguments to pass to plot, just so we don't keep repeating them
        kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
        ax.plot((-d, +d), (-d, +d), **kwargs)        # top-left diagonal
        ax2 = ax_bottom
        kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
        ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
        #remove one of the legend
        ax_bottom.legend_.remove()
        ax_top.set_xlabel('')
        plt.xlabel('')
        ax_top.tick_params(bottom = False)
        ax_top.set_ylabel('')
        ax_bottom.set_ylabel('')
        f.text(0.04, 0.5, 'Residence time (s)', ha='center', va='center', rotation='vertical')
        ax_top.legend(title = '')
        f.savefig(f"{plot_export}/Violin_plot_splitaxis.svg", dpi = 600)
        plt.show()

plot_violin(final, 'split')

############# Code to generate and collate summary statistics of dwell times
mean = final.groupby(['treatment', 'transition_type']).mean()
sem =  final.groupby(['treatment', 'transition_type']).sem()
final_drop = final.drop('transition_name', axis = 1)
N = final_drop.groupby(['treatment', 'transition_type']).count()
collated = pd.concat([mean,sem, N], axis = 1)
collated.drop([col for col in collated.columns.tolist() if 'y_axis_log10' in col], axis = 1, inplace = True)
collated.columns = ['mean_residence_time', 'sem', 'n']
collated.reset_index(inplace = True)
collated.to_csv(f"{output_folder}/summary.csv", index = False)
collated_filt = collated[(collated['transition_type'] == f'< {FRET_thresh} to > {FRET_thresh}')|(collated['transition_type'] == f'> {FRET_thresh} to < {FRET_thresh}')]
final_filt = final[(final['transition_type'] == f'< {FRET_thresh} to > {FRET_thresh}')|(final['transition_type'] == f'> {FRET_thresh} to < {FRET_thresh}')]


def plot_bar_with_sem(df, summary_df, y_axis = 'y_axis', palette = 'mako', order = order):
    ###### prepare order of datasets so that sem are correctly mapped onto the correct dataset
    list_to_order = list(np.arange(0, len(order), 1))
    dict_to_order = dict(zip(order, list_to_order))
    summary_df['plot_order'] = summary_df['treatment'].map(dict_to_order)
    collated_sorted = summary_df.sort_values(['plot_order', 'transition_type'])
    sorted_df = df.sort_values(['treatment', 'transition_type'])
    ###### now plot the figure
    fig, ax = plt.subplots()
    sns.set(style = 'ticks', font_scale = 1)
    sns.barplot(x='transition_name', y=y_axis, data=sorted_df, hue= 'treatment', palette=palette, ci =None, hue_order=order)
    x_coords = [p.get_x() + 0.5 * p.get_width() for p in ax.patches]
    y_coords = [p.get_height() for p in ax.patches]
    ax.errorbar(x=x_coords, y=y_coords, yerr=collated_sorted["sem"], fmt="none",  elinewidth = 2, capsize = 4, color = 'black')
    plt.ylabel('Residence time (s)')
    plt.xlabel('')
    plt.legend(title = '')
    fig.savefig(f'{plot_export}/mean_residence_withSEM_{y_axis}.svg', dpi = 600)
    plt.show()



def plot_residence_time_of_class(df, binwidth, transition, plot_type = 'KDE', log = False):
    if plot_type == 'KDE':
        for transition, dfs in df.groupby('transition_type'):
            filt_trans = dfs[dfs['transition_type'] ==  transition]
            sns.kdeplot(data = filt_trans, x = 'y_axis', hue = 'treatment', fill = False,  log_scale = True, common_norm=False, palette = 'mako')
            plt.title(transition)
            plt.show()
    if plot_type == 'individual':
        for treatment, dfs in df.groupby('treatment'):
            fig, axes = plt.subplots(4, 1, figsize = (8, 18), sharex=True)
            axes = axes.flatten()
            for i, transition in enumerate(list(dfs['transition_type'].unique())):
                fig = sns.histplot(data = dfs[dfs['transition_type']==transition], x = 'y_axis', binwidth = binwidth, kde = True, stat = 'density', log_scale = log, ax = axes[i])
                axes[i].set_xlabel("Residence time before transition to 'bound' state (s)")
                axes[i].set_title(f'{treatment} and {transition}')
                plt.xlim(0, 50)
            plt.savefig(f'{plot_export}/residence_time_histogram_{plot_type}.svg', dpi = 600)
            plt.show()
    if plot_type == 'collated':
        fig, axes = plt.subplots(5, 1, sharex = True)
        axes = axes.flatten()
        for i, treatment in enumerate(list(df['treatment'].unique())):
            dfs = df[df['treatment']==treatment]
            df2 = dfs[dfs['transition_type']==transition]
            fig = sns.histplot(
                data = df2, 
                hue = 'treatment', 
                x = 'y_axis',
                binwidth = binwidth, 
                # kde = True, 
                stat = 'density', 
                log_scale = log, 
                ax = axes[i], 
                common_norm = False, 
                fill = False, 
                palette = 'mako')
            axes[i].set_xlabel("Residence time before transition to 'bound' state (s)")
            axes[i].set_title(f'{treatment}', loc = 'center')
            plt.xlim(0, 200)
        plt.savefig(f'{plot_export}/residence_time_histogram_{plot_type}.svg', dpi = 600)
        plt.show()
    if plot_type == 'cum_dwell':
        fig, axes = plt.subplots(5, 1, sharex = True)
        axes = axes.flatten()
        for i, treatment in enumerate(list(df['treatment'].unique())):
            dfs = df[df['treatment']==treatment]
            df2 = dfs[dfs['transition_type']==transition]
            fig = sns.histplot(
                data = df2, 
                hue = 'treatment', 
                x = 'CumulativeTime(s)',
                binwidth = binwidth, 
                kde = True, 
                stat = 'density', 
                log_scale = log, 
                ax = axes[i], 
                common_norm = False, 
                fill = False, 
                palette = 'mako')
            axes[i].set_xlabel("Residence time before transition to 'bound' state (s)")
            axes[i].set_title(f'{treatment}')
            axes[i].legend('')
            plt.xlim(0, 200)
        plt.savefig(f'{plot_export}/residence_time_histogram_{plot_type}.svg', dpi = 600)
        plt.show()
    plt.show



###################
################### plot regular histogram data
##################


plot_bar_with_sem(final, collated, 'y_axis','BuPu', order)
plot_bar_with_sem(final_filt, collated_filt, 'y_axis', 'BuPu', order)
plot_residence_time_of_class(final, 20, f'< {FRET_thresh} to > {FRET_thresh}','collated', False)

#################
################# Prepare cumulative dwell time data for plotting
#################


dict_for_label_cum_dwell = {
'low-low':'$T_{low-low}$', 
'low-high':'$T_{low-high}$',
'high-high':'$T_{high-high}$',
'high-low':'$T_{high-low}$'
}

cumulative_dwell = cumulative_dwell[cumulative_dwell['treatment_name'].isin(order)]
cumulative_dwell['transition_name'] = cumulative_dwell['transition_type'].map(dict_for_label_cum_dwell)
mean = cumulative_dwell.groupby(['treatment_name', 'transition_type'])['CumulativeTime(s)'].mean()
sem =  cumulative_dwell.groupby(['treatment_name', 'transition_type'])['CumulativeTime(s)'].sem()
cumulative_dwell_drop = cumulative_dwell.drop('transition_name', axis = 1)
N = cumulative_dwell_drop.groupby(['treatment_name', 'transition_type'])['CumulativeTime(s)'].count()
col_cum_dwell = pd.concat([mean,sem, N], axis = 1)
col_cum_dwell.drop([col for col in col_cum_dwell.columns.tolist() if 'y_axis_log10' in col], axis = 1, inplace = True)
col_cum_dwell.columns = ['mean_residence_time', 'sem', 'n']
col_cum_dwell.reset_index(inplace = True)
cumulative_dwell.columns = ['molecule', 'FRET_before', 'FRET_after', 'Time', 'number_of_frames', 'treatment', 'transition_type', 'shift', 'is_in_sequence', 'CumulativeTime', 'CumulativeTime(s)', 'transition_name']
col_cum_dwell.columns = ['treatment', 'transition_type', 'mean_residence_time', 'sem', 'n']
col_cum_dwell.to_csv(f"{output_folder}/summary_cum_dwell.csv", index = False)
col_cum_dwell_filt = col_cum_dwell[(col_cum_dwell['transition_type'] == 'low-high')|(col_cum_dwell['transition_type'] == 'high-low')]
cumulative_dwell_filt = cumulative_dwell[(cumulative_dwell['transition_type'] == 'low-high')|(cumulative_dwell['transition_type'] == 'high-low')]

###################
################### plot cumulative residence time data
##################

plot_bar_with_sem(cumulative_dwell, col_cum_dwell, 'CumulativeTime(s)', 'BuPu', order)
plot_bar_with_sem(cumulative_dwell_filt, col_cum_dwell_filt, 'CumulativeTime(s)','BuPu', order)
plot_residence_time_of_class(cumulative_dwell_filt, 20, 'low-high','cum_dwell', False)

