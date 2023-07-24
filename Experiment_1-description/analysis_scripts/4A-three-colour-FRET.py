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

output_folder = "Experiment_1-RNA09-AF488/raw_data/"  ### Change for each experiment
plot_export = 'Experiment_1-RNA09-AF488/python_results/Traces-PPR-FRET/'
if not os.path.exists(plot_export):
    os.makedirs(plot_export)

input_folders = [
'Z:/Chaperone_subgroup/NickM/PPR/TIRF/230620_RNA09-AF488/Chan1_RNA09-AF488_20nM_steady/Channel532 - andor,488 - 200ms_Seq0000/Chan1_RNA09-AF488_20nM_steady_BC_all/', 
'Z:/Chaperone_subgroup/NickM/PPR/TIRF/230620_RNA09-AF488/Chan2_RNA09-AF488_20nM_duplexto-DNA02/Channel532 - andor,488 - 200ms_Seq0000/RNA09-DNA02-duplex_pprfret/', 
'Z:/Chaperone_subgroup/NickM/PPR/TIRF/230620_RNA09-AF488/Chan3_RNA01-AF488_20nM_steady/Channel532 - andor,488 - 200ms_Seq0000/RNA01-AF488_20nM_pprfret/', 

'Z:/Chaperone_subgroup/NickM/PPR/TIRF/230620_RNA09-AF488/Chan3_RNA01-AF488_20nM_steady/Channel532 - andor,488 - 200ms_Seq0000/RNA01-AF488_combined-with-230214_pprfret/', 
]

def move_folders(input_folders, filetype, output_folder):
    for folder in input_folders:
        new_folder = folder.split("/")[-2]
        if not os.path.exists(f"{output_folder}{new_folder}/"):
            os.makedirs(f"{output_folder}{new_folder}/")
        filelist = [filename for filename in os.listdir(folder) if filetype in filename] # create list of files
        for filename in filelist:
            shutil.copyfile(f"{folder}{filename}", f"{output_folder}{new_folder}/{filename}")        # if need to change filenames could have a function defined above that changes it
move_folders(input_folders, ".txt", output_folder)

dict_to_concat = {
    'RNA09-AF488':'Experiment_1-RNA09-AF488/raw_data/Chan1_RNA09-AF488_20nM_steady_BC_all',
    # 'RNA09-duplex':'Experiment_1-RNA09-AF488/raw_data/RNA09-DNA02-duplex_pprfret',
    # 'RNA01-AF488':'Experiment_1-RNA09-AF488/raw_data/RNA01-AF488_20nM_pprfret',
    'RNA01-AF488_col':'Experiment_1-RNA09-AF488/raw_data/RNA01-AF488_combined-with-230214_pprfret',

 }

font = {'family' : 'arial',
'weight' : 'normal'}
plt.rc('font', **font)
plt.rcParams['svg.fonttype'] = 'none'

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


# compiled_df[(compiled_df['molecule number']=='13')&(compiled_df['treatment']=='RNA01-DNA01-duplex-AF488')]
#   #### .rename(columns = {1:"test", 3:"test2"}) ## can rename individually if needed

for (treatment, molecule), df in compiled_df.groupby(['treatment','cumulative_molecule']):
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


def plot_FRET_multiple(df, top_FRET, bottom_FRET):
    data_hist = df[(df[top_FRET] > -0.2) & (df[top_FRET] < 1.2) & (df[bottom_FRET] > -0.2) & (df[bottom_FRET] < 1.2)]
    fig, axes = plt.subplots(2, 1, sharex =  True)
    sns.set_style("whitegrid",{'figure.figsize':(8,5.5), 'grid.linestyle':'--', 'font_scale':1.5} )
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


plot_FRET_multiple(compiled_df, 'FRET Cy3 to AF647', 'FRET Cy3 to AF647')
plot_FRET_multiple(FRET_bound, 'Cy3 FRET cascade', 'AF647 FRET cascade')

RNA01_df = compiled_df[compiled_df['treatment'] == 'RNA01-AF488_col']
RNA09_df = compiled_df[compiled_df['treatment'] == 'RNA09-AF488']

def plot_hist_type(df, kind = 'kde'):
    plot_hist, ax = plt.subplots()
    sns.set_style("ticks",{'font_scale':1} )
    if kind == 'kde':
        sns.kdeplot(
            data = df, 
            palette = 'greys', 
            x = "FRET Cy3 to AF647",
            hue="treatment",
            common_norm=False, 
            fill = True, 
            linewidth = 1.5, 
            alpha = 0.25)
    if kind == 'bar':
        sns.histplot(
            data = df, 
            x = "FRET Cy3 to AF647",
            common_norm=False, 
            stat = 'density',
            # hue = 'treatment',
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
    plot_hist.savefig(f'{plot_export}/Histogram_{kind}.svg', dpi = 600)
    plt.show()

plot_hist_type(compiled_df, 'bar')
plot_hist_type(RNA01_df, 'bar')
plot_hist_type(RNA09_df, 'bar')




palette_dict_treatment = {
    'RNA01-AF488':{'RNA-bound':'black', 'Unbound':'grey'}, 
    'RNA09-AF488':{'RNA-bound':'royalblue', 'Unbound':'skyblue'}, 
    'RNA09-duplex':{'RNA-bound':'purple', 'Unbound':'orchid'}, 
    'RNA01-AF488_col':{'RNA-bound':'black', 'Unbound':'grey'}, 


    }
palette = {True:'rebeccapurple', False:'orange'}

def plot_intensities1(df, intensity_type, palette):
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

fig, axes = plot_intensities1(compiled_df, 'AF647 at 488', palette_dict_treatment)


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
        # handles, labels = plt.gca().get_legend_handles_labels()
        # by_label = dict(zip(labels, handles))
        # axes[i].legend(['Unbound', 'Bound'], [by_label[key] for key in ['Unbound', 'RNA-bound']])
    axes[2].set_xlabel('Total fluorescence (a.u.)')
    plt.tight_layout()
    plt.xlim(-5000, 20000)
    fig.savefig(f'{plot_export}/intensity_bound_{treatment}.svg', dpi = 600)
    plt.show()
    return fig, axes


for treatment, df in compiled_df.groupby('treatment'):
    fig, axes = plot_each_intensity(df, treatment)



