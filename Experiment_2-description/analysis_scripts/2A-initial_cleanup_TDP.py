from scipy.stats import norm
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import pandas as pd
import numpy as np

output_folder = "Experiment_X-description/python_results"  ## change for each result

#### the first item in the tuple will be the name that goes into the graph legend
data_paths = {
    "treatment_label_1":"data_directory_1",
    "treatment_label_2":"data_directory_1",
    "treatment_label_3":"data_directory_1",
    "treatment_label_4":"data_directory_1",
    "treatment_label_5":'data_directory_1'
}

from Utilities.Data_analysis import remove_outliers, filter_TDP, file_reader, count_filt_mol

compiled_df = []
for data_name, data_path in data_paths.items():
    imported_TDPdata = file_reader(data_path, 'TDP')
    cleaned_TDPdata = remove_outliers(imported_TDPdata, 'TDP')
    cleaned_TDPdata["treatment name"] = data_name
    compiled_df.append(cleaned_TDPdata)
compiled_df = pd.concat(compiled_df)   #### .rename(columns = {1:"test", 3:"test2"}) ## can rename individually if needed
compiled_df.columns = ['Molecule', 'FRET before transition', 'FRET after transition', 'Time', "treatment_name"]


filtered_data = filter_TDP(compiled_df, 0.3)  ##### number refers to the FRET threshold to filter data. Will only include molecules that go below this set threshold. Will default to 0.3 
filtered_data.to_csv(f'{output_folder}/TDP_cleaned_filt.csv', index = False)
############### Save data to TDP_cleaned
compiled_df.to_csv(f'{output_folder}/TDP_cleaned.csv', index = False)

############### calculate mol count
FRET_value = 0.6
molcount = count_filt_mol(compiled_df, FRET_value, data_paths, 1)  ## last number refers to the index of the treatment you want to subtract from others
molcount.to_csv(f'{output_folder}/mol_below_{FRET_value}.csv', index = None)


