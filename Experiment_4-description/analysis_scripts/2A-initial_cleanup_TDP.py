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

def file_reader(input_folder):
    filename = input_folder
    A = pd.read_table(filename, header = None)
    A.columns = ['Molecule', 'FRET before transition', 'FRET after transition', 'Time']
    return A

def remove_outliers(compiled_TDP):
    outliers = compiled_TDP[(compiled_TDP["FRET before transition"] < -0.5)|(compiled_TDP["FRET before transition"] > 1.5)|(compiled_TDP["FRET after transition"] < -0.5) | (compiled_TDP["FRET after transition"] > 1.5)].index
    compiled_TDP.drop(outliers, inplace = True)
    return compiled_TDP

def filter_TDP(data_frame):
    filtered_mol = []
    for treatment, df in data_frame.groupby("treatment_name"):
        mol_list = df[(df["FRET before transition"] <= 0.3)|(df["FRET after transition"] <= 0.3)].Molecule.unique().tolist()
        filtered = df[df["Molecule"].isin(mol_list)]
        filtered_mol.append(filtered)
    filtered_mol = pd.concat(filtered_mol)
    return filtered_mol


compiled_df = []
for data_name, data_path in data_paths.items():
    imported_TDPdata = file_reader(data_path)
    cleaned_TDPdata = remove_outliers(imported_TDPdata)
    cleaned_TDPdata["treatment name"] = data_name
    compiled_df.append(cleaned_TDPdata)
compiled_df = pd.concat(compiled_df)   #### .rename(columns = {1:"test", 3:"test2"}) ## can rename individually if needed
compiled_df.columns = ['Molecule', 'FRET before transition', 'FRET after transition', 'Time', "treatment_name"]


filtered_data = filter_TDP(compiled_df)
filtered_data.to_csv(f'{output_folder}/TDP_cleaned_filt.csv', index = False)
############### Save data to TDP_cleaned
compiled_df.to_csv(f'{output_folder}/TDP_cleaned.csv', index = False)

