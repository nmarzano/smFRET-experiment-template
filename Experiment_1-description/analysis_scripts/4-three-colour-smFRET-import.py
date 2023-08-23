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

output_folder = "Experiment_2-description/raw_data/"  ### Change for each experiment

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

