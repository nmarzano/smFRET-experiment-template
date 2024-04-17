from smfret.src.plotting_scripts import histogram_plots as pa
from smfret.src.plotting_scripts import timelapse_plots as pb
from smfret.src.plotting_scripts import traces_plot as pc
from smfret.src.plotting_scripts import heatmap_liveflow_plots as pe
from smfret.src.plotting_scripts import heatmap_liveflow_plots as pf
from smfret.src.plotting_scripts import synchronised_transition_plots as pg
from smfret.src.plotting_scripts import gaussian_plots as ph

import glob as glob
import os as os


if __name__ == "__main__":

    output_folder = "Experiment_1-description/python_results"

    data_paths = {
        "6min":("6min", "treatment/"),
        "12min":("12min", "treatment/"),
        "18min":("18min", "treatment/"),
        "24min":("24min", "treatment/"),
        "30min":("30min", "treatment/"),
        "36min":("36min", "treatment/"),
        "42min":("42min", "treatment/"),
        "48min":("48min", "treatment/"),
        "54min":("54min", "treatment/"),
        "60min":("60min", "treatment/"),
    }
    # ---------------------- plot histograms and define color scheme --------------------------------------------------------


    colors = {
        "Native": "black", 
        "Spontaneous":"darkorange", 
        "treatment_3":"skyblue", 
        "treatment_4":"royalblue", 
        "treatment_5":"darkviolet"
    }

    compiled_df, filt_dfs = pa.master_histogram_func(data_paths, output_folder=output_folder, thresh=0.3)

    # -------------------- combine timelapse data from different treatments and scripts -----------------------------------------------

    data_from_exp = {
    'treatment1':'Experiment_1-description/python_results/mean.csv',
    'treatment2':'Experiment_2-description/python_results/mean.csv',
    'treatment3':'Experiment_3-description/python_results/mean.csv',
    'treatment4':'Experiment_4-description/python_results/mean.csv',
    }

    timelapse_colors = {
    'treatment1':'black',    
    'treatment2':'purple', 
    'treatment3':'skyblue',
    'treatment4':'royalblue'
    }

    pb.master_timelapse_func(data_from_exp, 
                             thresh=0.3, 
                             xlim_min=6, 
                             xlim_max=60, 
                             output_folder=output_folder, 
                             palette=timelapse_colors, 
                             data_type='Mean')  #### can either be raw values 'Mean' or normalised to 100% 'normalised'.


    # ------------------- plot all traces with defined length or specific molecules of interest -------------------------------------
    
    data_paths_traces = {
        "example1":"trace1.dat",
        "example2":"trace2.dat",
    }

    pc.master_plot_individual_trace(data_paths_traces, input_folder=output_folder, exposure=0.2)
    pc.master_plot_traces_func(input_folder=output_folder, exposure=0.2, min_trace_length=180)

    # ------------------ process and plot data for live flow-in experiments -------------------------
    data_paths_liveflow = {
        'treatment':('treatment_name', 'Experiment_2-description/raw_data/treatment_to_plot'),
        'treatment':('treatment_name', 'Experiment_2-description/raw_data/treatment_to_plot'),
    }

    transition_type = 'low_to_high' ### or high_to_low

    pe.master_heatmap_processing_func(data_paths_liveflow, 
                                input_folder=output_folder, 
                                exposure=0.200, 
                                FRET_thresh=0.3, 
                                transition_type=transition_type, 
                                time_thresh=15, 
                                injection_time=15)

    

    pf.master_plot_flowin_func(input_folder=output_folder, 
                            transition_type=transition_type, 
                            gridsize=100, 
                            binshex=80)



    # --------------------------------- synchronised transitions ------------------------------------------
    order = list(data_paths.keys())
    list_to_drop = ['']
    filtered_list = [i for i in order if i not in list_to_drop]

    percent_trans_meet_criteria_df, calculated_transitions_df, consecutive_from_dnak_release, nonconsecutive_from_dnak_release = pg.master_plot_synchronised_transitions(order=order, 
                                        output_folder=output_folder, 
                                        exposure=0.2, 
                                        frames_to_plot=50, 
                                        FRET_before=0.3, 
                                        FRET_after=0.3, 
                                        datatype="Proportion", 
                                        filt=False, 
                                        palette='BuPu')#### datatype could be ratio



    # --------------------------- gaussian fitting and plottting proportion of populations over time -------------------------------

    
    timepoint = {'6min':0, '12min':6,'18min':12,'24min':18,'30min':24,'36min':30}

    collated_df, gauss_saveloc = ph.fit_gauss_master(input_folder=output_folder, gauss_num='three', timepoint=None)
    test, collated_df, gauss_saveloc = ph.fit_gauss_master(input_folder=output_folder, gauss_num='three', timepoint=timepoint)

    # Combined the gaussian fit proportion data from multiple experiments and combine here to plot

    data_from_exp = {
    'KJE':'Experiment_1-KJE-timelapse/python_results/GaussianFits/collated_populations.csv',
    'KJE-lowG':'Experiment_5-low-G-timelapse/python_results/GaussianFits/collated_populations.csv',
    'KJE-highG':'Experiment_2-KJEG-timelapse/python_results/GaussianFits/collated_populations.csv',
    }

    final = ph.plot_gauss_timelapse(data_from_exp, gauss_saveloc)



