import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import pandas as pd
import matplotlib
from scipy import optimize, signal
import numpy as np
from scipy import optimize
import scipy.integrate as integrate
from lmfit import models
import Utilities.Data_analysis as uda

output_folder = "Experiment_5-IDS1_KJE-timelapse/python_results"  ### modify for each experiment
thresh = 0.2


data_from_exp = {
'IDS1_KJE':f'{output_folder}/mean.csv',
'NDS1_KJE':'Experiment_4-NDS1_KJE-only_timelapse/python_results/mean.csv',
'NDS1_KJEG':'Experiment_3-NDS1_col_timelapse/python_results/mean.csv',
'IDS1_KJEG':'Experiment_2-IDS1_col_timelapse/python_results/mean.csv',
}

colors_plot = {
'NDS1_KJE':'black', 
'NDS1_KJEG':'purple', 
'IDS1_KJEG':'skyblue',
'IDS1_KJE':'royalblue'
}


###### Equation for exponential, can always define new equation and call in the 
###### plot data col function

def guess_exponential(x,A,B):
    y = A*np.exp(-B*x)
    return y


########## Generic function used to fit curve to plot. The fit should be easily changed by
########## calling a different equation into 'fit_type).

def fit_curve_to_plot(df, fit_type, x, y, data = 'Mean'):
    est_x = df['timepoint'].to_numpy()
    est_y = df[f'{data}'].to_numpy()
    # actual code for fitting
    paramaters, covariance = optimize.curve_fit(fit_type, est_x, est_y)
    fit_A = paramaters[0]
    fit_B = paramaters[1]
    print(fit_A)
    print(fit_B)
    #plotting data from fit
    fit_x_values = np.linspace(x,y,1000)
    fit_y = fit_type(fit_x_values, fit_A, fit_B)
    return fit_x_values, fit_y


##############
############## Code to combine timepoint data from multiple mutants/conditions and plot in a single graph
##############

test = []
for data_name, data_path in data_from_exp.items():
    data = uda.file_reader(data_path, 'other')
    data['treatment'] = data_name
    test.append(data)
final = pd.concat(test).reset_index()
final = final.iloc[:, 2:]


def plot_data_col(df, fit_type, xlim, ylim, data = 'Mean'):
    fig, ax = plt.subplots()
    for treatment, dfs in df.groupby('treatment'):
        if data == 'Mean':
            fitx, fity = fit_curve_to_plot(dfs, fit_type, xlim, ylim)
        else:
            fitx, fity = fit_curve_to_plot(dfs, fit_type, xlim, ylim, 'normalised')
        sns.scatterplot(data = dfs, x = 'timepoint', y = data, hue='treatment', palette=colors_plot)
        plt.errorbar(dfs['timepoint'],dfs[f'{data}'],dfs['Std_error'], fmt='none', capsize=3, ecolor='black')
        plt.plot(fitx, fity,'k', color=colors_plot[f'{treatment}'])
    plt.ylabel(f'Proportion of time spent < {thresh} FRET (mean)')
    plt.xlabel('Time (min)')
    plt.legend(title = '')
    plt.savefig(f'{output_folder}/timelapse_proportion_mol_below_thresh.svg', dpi = 600)
    plt.show()
    return 

plot_data_col(final, guess_exponential, 6, 36, 'normalised')
