import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import pandas as pd
import matplotlib

output_folder = "Experiment_X-description/python_results"  ### modify for each experiment

#### the first item in the tuple will be the name that goes into the graph legend
data_paths = {
    "treatment_label_1":("figure_label", "data_directory_1/"),
    "treatment_label_2":("figure_label", "data_directory_2/"),
    "treatment_label_3":("figure_label", "data_directory_3/"),
    "treatment_label_4":("figure_label", "data_directory_4/"),
    "treatment_label_5":("figure_label", "data_directory_5/")
}

from Utilities.Data_analysis import remove_outliers, file_reader

compiled_df = []
for data_name, (label, data_path) in data_paths.items():
    imported_data = file_reader(data_path, 'hist')
    cleaned_raw = remove_outliers(imported_data, 'hist')    #### add "idealized" after imported_data to get idealized histograms
    cleaned_raw["treatment_name"] = data_name
    compiled_df.append(cleaned_raw)
compiled_df = pd.concat(compiled_df)   #### .rename(columns = {1:"test", 3:"test2"}) ## can rename individually if needed
compiled_df.columns = ["frames", "donor", "acceptor", "FRET", "idealized FRET", "treatment_name"]
compiled_df.to_csv(f'{output_folder}/Raw_FRET_histogram_data.csv', index = False)

################## Code to assemble relevant data into an array and then plot using a single line of code
colors = {
    "Native": "black", 
    "Spontaneous":"darkorange", 
    "treatment_3":"skyblue", 
    "treatment_4":"royalblue", 
    "treatment_5":"darkviolet"
}

labels = {data_name:label for data_name, (label, data_path) in data_paths.items()}

order = ["KJE_late","KJE_early", "KJ_ATP", "Spontaneous", "Native"]


font = {'weight' : 'normal',
'size'   : 12 }
matplotlib.rcParams['font.sans-serif'] = "Arial"
matplotlib.rcParams['font.family'] = "sans-serif"
matplotlib.rc('font', **font)
plt.rcParams['svg.fonttype'] = 'none'
sns.set(style="darkgrid", font_scale = 1.5, rc={'figure.figsize':(8,5.5)})
plot_hist = sns.kdeplot(
    data = compiled_df, 
    palette = colors, 
    x = "FRET",
    hue="treatment_name",
    hue_order = order,
    common_norm=False, 
    fill = True, 
    linewidth = 1.5, 
    alpha = 0.25)
plt.xlim(0, 1, 10)
plt.xlabel("FRET")
plt.legend([labels[treatment] for treatment in reversed(order)],loc='upper left', fontsize = 15)

plot_hist.figure.savefig(f'{output_folder}/Histogram.svg', dpi = 600)
plt.show()


