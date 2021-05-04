from scipy.stats import norm
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import pandas as pd
import numpy as np
import glob


filename = 'Experiment_X-description/python_results/TDP_cleaned.csv'
output_folder = 'Experiment_X-description/python_results'


FRET_thresh = 0.5 #### FRET value at which to filter data above or below
fps = 5  ### frames per second
thresh = 2 ### should be 10x expsoure if using NL515 smoothing on MASH FRET
headers = ["< 0.5 to < 0.5", "< 0.5 to > 0.5", "> 0.5 to > 0.5", "> 0.5 to < 0.5"]
headers_withsum =  ["< 0.5 to < 0.5", "< 0.5 to > 0.5", "> 0.5 to > 0.5", "> 0.5 to < 0.5", "sum", "sample"]
TDP_data = pd.read_csv(filename, header = "infer")

#######################################################################################################################
####################################  Code to remove the first dwell time from each molecule - does this by making an empty list and then a "for" loop
####################################  Loop structured so that the dataframe "A" is grouped according to the molecule number in the column "Molecule"
####################################  Loop will go through each molecule in df (i.e. dataframe "A") and will subset all indexes except for the first one
####################################  (i.e iloc[1:]) and then append that data (now without the first dwell time) into the list made at the beginning. 
####################################  Data is then concatenated into a dataframe called filtered (overwrites the list)
def cleanup_dwell(data, first_dwell = "delete"):
    if first_dwell == "delete":
        filtered = []
        for molecule, df in data.groupby("Molecule"):
            filtered.append(df.iloc[1:])
        filtered = pd.concat(filtered)  #####filtered = pd.concat([df.iloc[1:] for molecule, df in A.groupby("Molecule")]) ##code here is the same as the for loop but in a list comprehension format
        filtered["Time (s)"] = filtered["Time"]/fps
        filtered = filtered[filtered["Time (s)"] >= thresh]
        return filtered
    if first_dwell == "keep":
        data["Time (s)"] = data["Time"]/fps
        data = data[data["Time (s)"] >= thresh]
        return data

################################### Code to start filtering FRET before and after values and extract the dwell times and then concatenate each result into
################################### a new data frame, df_col

def filter_dwell(df):
    filtered_lowtohigh = df[(df["FRET before transition"] < FRET_thresh) & (df["FRET after transition"] > FRET_thresh)].copy()
    filtered_lowtolow = df[(df["FRET before transition"] < FRET_thresh) & (df["FRET after transition"] < FRET_thresh)].copy()
    filtered_hightolow = df[(df["FRET before transition"] > FRET_thresh) & (df["FRET after transition"] < FRET_thresh)].copy()
    filtered_hightohigh = df[(df["FRET before transition"] > FRET_thresh) & (df["FRET after transition"] > FRET_thresh)].copy()
    dataf = [filtered_lowtolow["Time (s)"], filtered_lowtohigh["Time (s)"], filtered_hightohigh["Time (s)"], filtered_hightolow["Time (s)"]]
    df_col = pd.concat(dataf, axis = 1, keys = headers)
    df_col = df_col.apply(lambda x:pd.Series(x.dropna().values))  ## removes NaN values from each column in df_col
    return df_col


################################### Calculate the percentage probability of each transition occuring
def transition_frequency(filt):
    count_df_col = pd.DataFrame(filt.count(axis = 0)).transpose()
    count_df_col["sum"] = count_df_col.sum(axis = 1)
    dwell_frequency = pd.DataFrame([(count_df_col[column]/count_df_col["sum"])*100 for column in count_df_col]).transpose()
    print(dwell_frequency)
    return dwell_frequency


################################### Calculate the mean of each transition type
def calculate_mean(filtered_data, treatment_name):
    mean_dwell = pd.DataFrame([filtered_data.iloc[0:].mean()])
    mean_dwell["sample"] = treatment_name
    return mean_dwell


for treatment_name, df in TDP_data.groupby("treatment_name"):
    initial_data = df[df["treatment_name"] == treatment_name]
    cleaned_data = cleanup_dwell(initial_data) ##### to keep the first dwell state, simply change code to "cleanup_dwell(initial_data, "keep")
    filtered_data = filter_dwell(cleaned_data)
    filtered_data.to_csv(f"{output_folder}/Dwell_times/Filtered_dwelltime_{treatment_name}.csv", index = False)
    mean_dwell = calculate_mean(filtered_data, treatment_name)
    mean_dwell.to_csv(f"{output_folder}/Mean_dwell/Filtered_meandwell_{treatment_name}.csv", index = False)
    dwell_frequency = transition_frequency(filtered_data)
    dwell_frequency["sample"] = treatment_name
    dwell_frequency.to_csv(f"{output_folder}/Dwell_frequency/Filtered_dwellfrequency_{treatment_name}.csv", index = False, header = None)
