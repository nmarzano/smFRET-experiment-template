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
compiled_df.columns = ["frames", "donor", "acceptor", "FRET", "idealized FRET", 'molecule_number', "treatment_name"]
compiled_df.to_csv(f'{output_folder}/Cleaned_FRET_histogram_data.csv', index = False)

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

###### FOR RIDGELINE PLOTS
###################################################
###################################################
data_paths_ridgeline = {
    "1":"label",
    '2':'label',
    "3":"label",
    "4":"label", 
    "5":"label", 
    "6":"label",
    "7":"label", 
    "8":"label", 
    "9":"label", 
    "10":"label", 
    "11":"label", 
    "12":"label",
}

pal = sns.color_palette(palette='mako', n_colors=12)
matplotlib.rc('font', **font)
plt.rcParams['svg.fonttype'] = 'none'
g = sns.FacetGrid(compiled_df, row='treatment_name', hue='treatment_name', aspect=10, height=5, palette=pal)
# then we add the densities kdeplots for each condition
g.map(sns.kdeplot, 'FRET',
      bw_adjust=1, clip_on=False,
      fill=True, alpha=1, linewidth=10)
# here we add a white line that represents the contour of each kdeplot
g.map(sns.kdeplot, 'FRET', 
      bw_adjust=1, clip_on=False, 
      color="white", linewidth=6)
# here we add a horizontal line for each plot
g.map(plt.axhline, y=0,
      lw=2, clip_on=False)
# we loop over the FacetGrid figure axes (g.axes.flat) and add the condition as text with the right color
# notice how ax.lines[-1].get_color() enables you to access the last line's color in each matplotlib.Axes
for i, ax in enumerate(g.axes.flat):
    ax.text(0, .5, data_paths_ridgeline[f'{i+1}'],
            fontweight='bold', fontsize=100,
            color=ax.lines[-1].get_color())
    ax.set_facecolor((0, 0, 0, 0))  ### removes background so KDE plots are not obstructed when stacked
# we use matplotlib.Figure.subplots_adjust() function to get the subplots to overlap
g.fig.subplots_adjust(hspace=-0.7)
# eventually we remove axes titles, yticks and spines
g.set_titles("")
g.set(yticks=([]))
g.despine(bottom=True, left=True)
plt.setp(ax.get_xticklabels(), fontsize=100, fontweight='bold')
plt.xlabel('FRET', fontweight='bold', fontsize=100)
plt.xlim(-.4, 1.2)
g.savefig(f'{output_folder}/Histogram-ridgeline.svg', dpi = 600)
plt.show()


def plot_hist(df, order):
    sns.set(style="darkgrid", font_scale = 1.5, rc={'figure.figsize':(8,5.5)})
    plot_hist = sns.kdeplot(
        data = df, 
        palette = 'mako', 
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

plot_hist(compiled_df, order)



def fit_gauss_dif_constrained_allpeaks(df, treatment, mu_1, sigma_1, amplitude_1, gamma_1, mu_2, sigma_2, amplitude_2, mu_3, sigma_3, amplitude_3, gamma_3):
    filt_df = df[df['treatment_name'] == treatment]
    bins = np.arange(-0.21, 1.1, 0.025) 
    inds = np.digitize(filt_df['FRET'].astype(float), bins)
    xdata, ydata = np.unique(inds, return_counts=True)
    ydata = ydata[1:-1] #### trim off outside range bins at the end
    xdata = [np.mean(bins[x : x + 2]) for x in range(len(bins)- 1)]  ##### convert bin edges to bin centres, therefore end up with one less bin
    sns.lineplot(xdata, ydata)

    model_1 = models.SkewedGaussianModel(prefix='m1_')
    model_2 = models.GaussianModel(prefix='m2_')
    model_3 = models.SkewedGaussianModel(prefix='m3_')
    model = model_1 + model_2 + model_3 

    model_1.set_param_hint('m1_gamma', vary=False)
    model_1.set_param_hint('m1_sigma', vary=False)
    model_1.set_param_hint('m1_center', vary=False)

    model_2.set_param_hint('m2_sigma', vary=False)
    model_2.set_param_hint('m2_center', vary=False)
    model_3.set_param_hint('m3_gamma', vary=False)
    model_3.set_param_hint('m3_sigma', vary=False)


    params_1 = model_1.make_params(center = mu_1, sigma = sigma_1, amplitude = amplitude_1, gamma = gamma_1, min = 0)
    params_2 = model_2.make_params(center = mu_2, sigma = sigma_2, amplitude = amplitude_2, min = 0)
    params_3 = model_3.make_params(center = mu_3, sigma = sigma_3, amplitude = amplitude_3, gamma = gamma_3, min = 0)
    params = params_1.update(params_2)
    params = params.update(params_3)

    output = model.fit((ydata/np.max(ydata)), params, x=xdata)
    fig = output.plot(data_kws={'markersize': 3})

    paramaters = {name:output.params[name].value for name in output.params.keys()}
    fitx = np.arange(-0.2, 1.2, 0.025)

    fit1 = model_1.eval(x = fitx, center = paramaters['m1_center'], amplitude = abs(paramaters['m1_amplitude']), sigma = paramaters['m1_sigma'], gamma = paramaters['m1_gamma'])
    fit2 = model_2.eval(x = fitx, center = paramaters['m2_center'], amplitude = abs(paramaters['m2_amplitude']), sigma = paramaters['m2_sigma'], fwhm = paramaters['m2_fwhm'])
    fit3 = model_3.eval(x = fitx, center = paramaters['m3_center'], amplitude = abs(paramaters['m3_amplitude']), sigma = paramaters['m3_sigma'], gamma = paramaters['m3_gamma'])

    sns.lineplot(fitx, fit1)
    sns.lineplot(fitx, fit2)
    sns.lineplot(fitx, fit3)
    plt.show()

    # Calculate area under the curve for each gaussian
    aoc_m1 = paramaters['m1_amplitude']
    aoc_m2 = paramaters['m2_amplitude']
    aoc_m3 = paramaters['m3_amplitude']

    # aoc_m1 = (paramaters['m1_amplitude']*paramaters['m1_sigma'])/0.3989
    # aoc_m2 = (paramaters['m2_amplitude']*paramaters['m2_sigma'])/0.3989
    # aoc_m3 = (paramaters['m3_amplitude']*paramaters['m3_sigma'])/0.3989

    sum_aoc = aoc_m1 + aoc_m2 + aoc_m3 

    aoc_m1_percent_of_total = (aoc_m1/sum_aoc)*100
    aoc_m2_percent_of_total = (aoc_m2/sum_aoc)*100
    aoc_m3_percent_of_total = (aoc_m3/sum_aoc)*100
    list_of_gaus_proportion = [aoc_m1_percent_of_total, aoc_m2_percent_of_total, aoc_m3_percent_of_total]
    labels_of_gaus_proportion = ['m1', 'm2', 'm3']
    proportion_df = pd.DataFrame([labels_of_gaus_proportion, list_of_gaus_proportion])
    proportion_df.columns = proportion_df.iloc[0]
    proportion_df = proportion_df.drop(0)
    proportion_df.to_csv(f'{output_folder}/gaussian_proportions_for_{treatment}.csv')
    return proportion_df, paramaters


gaussian_medium_constrainall = fit_gauss_dif_constrained_allpeaks(compiled_df, 'medium', 0.02, .1, .5, 3, .63, .12, 1, .9, .22, .5, -2.7)
gaussian_col_high = fit_gauss_dif_constrained_allpeaks(compiled_df, 'high', 0.08, .08, .1, .1, .64, .10, .6, .9, .16, 1, -2.7)
gaussian_low_constrainall = fit_gauss_dif_constrained_allpeaks(compiled_df, 'low', 0.1, .1, .8, 2, .63, .11, .5, .95, .22, .4, -2.7)


def fit_gauss_constrain_native_center(df, treatment, mu_1, sigma_1, amplitude_1, gamma_1, mu_2, sigma_2, amplitude_2, mu_3, sigma_3, amplitude_3, gamma_3):
    filt_df = df[df['treatment_name'] == treatment]
    bins = np.arange(-0.21, 1.1, 0.025) 
    inds = np.digitize(filt_df['FRET'].astype(float), bins)
    xdata, ydata = np.unique(inds, return_counts=True)
    ydata = ydata[1:-1] #### trim off outside range bins at the end
    xdata = [np.mean(bins[x : x + 2]) for x in range(len(bins)- 1)]  ##### convert bin edges to bin centres, therefore end up with one less bin
    sns.lineplot(xdata, ydata)

    model_1 = models.SkewedGaussianModel(prefix='m1_')
    model_2 = models.GaussianModel(prefix='m2_')
    model_3 = models.SkewedGaussianModel(prefix='m3_')
    model = model_1 + model_2 + model_3 

    model_2.set_param_hint('m2_center', vary=False)
    model_3.set_param_hint('m3_gamma', vary=False)
    model_3.set_param_hint('m3_sigma', vary=False)

    params_1 = model_1.make_params(center = mu_1, sigma = sigma_1, amplitude = amplitude_1, gamma = gamma_1, min = 0)
    params_2 = model_2.make_params(center = mu_2, sigma = sigma_2, amplitude = amplitude_2, min = 0)
    params_3 = model_3.make_params(center = mu_3, sigma = sigma_3, amplitude = amplitude_3, gamma = gamma_3, min = 0)
    params = params_1.update(params_2)
    params = params.update(params_3)

    output = model.fit((ydata/np.max(ydata)), params, x=xdata)
    fig = output.plot(data_kws={'markersize': 3})

    paramaters = {name:output.params[name].value for name in output.params.keys()}
    fitx = np.arange(-0.2, 1.2, 0.025)

    fit1 = model_1.eval(x = fitx, center = paramaters['m1_center'], amplitude = abs(paramaters['m1_amplitude']), sigma = paramaters['m1_sigma'], gamma = paramaters['m1_gamma'])
    fit2 = model_2.eval(x = fitx, center = paramaters['m2_center'], amplitude = abs(paramaters['m2_amplitude']), sigma = paramaters['m2_sigma'], fwhm = paramaters['m2_fwhm'])
    fit3 = model_3.eval(x = fitx, center = paramaters['m3_center'], amplitude = abs(paramaters['m3_amplitude']), sigma = paramaters['m3_sigma'], gamma = paramaters['m3_gamma'])

    sns.lineplot(fitx, fit1)
    sns.lineplot(fitx, fit2)
    sns.lineplot(fitx, fit3)
    plt.show()

    # Calculate area under the curve for each gaussian
    aoc_m1 = paramaters['m1_amplitude']
    aoc_m2 = paramaters['m2_amplitude']
    aoc_m3 = paramaters['m3_amplitude']

    # aoc_m1 = (paramaters['m1_amplitude']*paramaters['m1_sigma'])/0.3989
    # aoc_m2 = (paramaters['m2_amplitude']*paramaters['m2_sigma'])/0.3989
    # aoc_m3 = (paramaters['m3_amplitude']*paramaters['m3_sigma'])/0.3989

    sum_aoc = aoc_m1 + aoc_m2 + aoc_m3 

    aoc_m1_percent_of_total = (aoc_m1/sum_aoc)*100
    aoc_m2_percent_of_total = (aoc_m2/sum_aoc)*100
    aoc_m3_percent_of_total = (aoc_m3/sum_aoc)*100
    list_of_gaus_proportion = [aoc_m1_percent_of_total, aoc_m2_percent_of_total, aoc_m3_percent_of_total]
    labels_of_gaus_proportion = ['m1', 'm2', 'm3']
    proportion_df = pd.DataFrame([labels_of_gaus_proportion, list_of_gaus_proportion])
    proportion_df.columns = proportion_df.iloc[0]
    proportion_df = proportion_df.drop(0)
    proportion_df.to_csv(f'{output_folder}/gaussian_proportions_for_{treatment}.csv')
    return proportion_df, paramaters

