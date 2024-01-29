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

output_folder = "Experiment_1-description/python_results"  ### modify for each experiment
thresh = 0.2
xlim_min = 6
xlim_max = 60

data_from_exp = {
'treatment1':f'{output_folder}/mean.csv',
'treatment2':'Experiment_2-description/python_results/mean.csv',
'treatment3':'Experiment_3-description/python_results/mean.csv',
'treatment4':'Experiment_4-description/python_results/mean.csv',
}

colors_plot = {
'treatment1':'black', 
'treatment2':'purple', 
'treatment3':'skyblue',
'treatment4':'royalblue'
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
    """Plots the proportion of time each molecule spends below a threshold (defined previously in 1A-plot-histogram) as the mean +- SE as a function of time. This function
    is predominantly designed to collate timelapse data from different experiments and present it within a single plot. It also fits a generic curve to the data, which has been 
    defined above in fit_curve_to_plot function. 

    Args:
        df (dataframe): collated dataset of all treatments to compare
        fit_type (func): the type of fit you wish to plot. If you want a different function, define elsewhere and call the function here to plot.
        xlim (float): minimum x-axis value used to define fit
        ylim (float): maximum x-axis value used to define fit
        data (str, optional): changes if you want to plot normalised data 'normalised' or raw data 'Mean'. Defaults to 'Mean'.
    """
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

plot_data_col(final, guess_exponential, xlim_min, xlim_max, 'normalised')  #### can either be raw values 'Mean' or normalised to 100% 'normalised'.
