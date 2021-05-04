from scipy.stats import norm
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import pandas as pd
import numpy as np
import glob

output_folder = "Experiment_X-description/python_results" #### NEED TO CHANGE EACH TIME TO THE PYTHON RESULTS LEVEL

headers_withsum =  ["< 0.5 to < 0.5", "< 0.5 to > 0.5", "> 0.5 to > 0.5", "> 0.5 to < 0.5", "sum", "sample"]
order = [4, 1, 2, 3, 0]  # organise data using its index 

def file_reader(input_folder): 
    filenames = glob.glob(input_folder + "/*.csv")

    dfs = []
    for filename in filenames:
        dfs.append(pd.read_csv(filename, header=None))
    test = pd.concat(dfs, ignore_index=True)
    test_dfs = pd.DataFrame(test)
    test_dfs.columns = headers_withsum
    return test_dfs

imported_df = file_reader(f'{output_folder}/Dwell_frequency') ##### NEED TO ALTER EACH TIME
imported_df.drop(columns = ["sum"], inplace = True)
imported_df = imported_df.reindex(order) ### need to do this to plot in order that you would like
imported_df.to_csv(f'{output_folder}/Transition_frequency.csv')


######################## code to prepare for plotting
total = [i+j+k+m for i,j,k, m in zip(imported_df['< 0.5 to < 0.5'], imported_df['< 0.5 to > 0.5'], imported_df['> 0.5 to > 0.5'], imported_df["> 0.5 to < 0.5"])]
lowtolow = [i / j * 100 for i,j in zip(imported_df['< 0.5 to < 0.5'], total)]
lowtohigh = [i / j * 100 for i,j in zip(imported_df['< 0.5 to > 0.5'], total)]
hightohigh = [i / j * 100 for i,j in zip(imported_df['> 0.5 to > 0.5'], total)]
hightolow = [i / j * 100 for i,j in zip(imported_df['> 0.5 to < 0.5'], total)]

################### plotting code
labels = imported_df["sample"].to_list()
barWidth = 0.85
plt.rcParams['svg.fonttype'] = 'none'
sns.set(style = "darkgrid")
# Create green Bars
plot1 = plt.figure()
plt.bar(labels, lowtolow, color='skyblue', edgecolor='white', width=barWidth, label = "< 0.5 to < 0.5" )
# Create orange Bars
plt.bar(labels, lowtohigh, bottom=lowtolow, color='royalblue', edgecolor='white', width=barWidth, label = "< 0.5 to > 0.5")
# Create blue Bars
plt.bar(labels, hightohigh, bottom=[i+j for i,j in zip(lowtolow, lowtohigh)], color='yellowgreen', edgecolor='white', width=barWidth, label = "> 0.5 to > 0.5")
# Create blue Bars
plt.bar(labels, hightolow, bottom=[i+j+k for i,j,k in zip(lowtolow, lowtohigh, hightohigh)], color='green', edgecolor='white', width=barWidth,label = "> 0.5 to < 0.5")
plt.legend(loc = "upper left", bbox_to_anchor = (1,1), ncol =1)
plt.ylabel("Transition probability (%)")
plot1.savefig(f'{output_folder}/Transition_frequency_plot.svg', dpi = 600)


