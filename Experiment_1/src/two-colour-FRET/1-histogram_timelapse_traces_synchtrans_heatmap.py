from smfret.src.plotting_scripts import histogram_plots as pa
from smfret.src.plotting_scripts import timelapse_plots as pb
from smfret.src.plotting_scripts import traces_plot as pc
from smfret.src.plotting_scripts import heatmap_liveflow_plots as pe
from smfret.src.plotting_scripts import heatmap_liveflow_plots as pf
from smfret.src.plotting_scripts import synchronised_transition_plots as pg
from smfret.src.plotting_scripts import gaussian_plots as ph
import glob as glob
import os as os
import pandas as pd

if __name__ == "__main__":

    output_folder = "Fig3_sensors_timelapse/python_results"
    FRET_thresh = 0.3
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
        "IDS1_KJE": "#c9ddf1", 
        "IDS1_KJEG":"#6999d1", 
        "IDS1_KJEGE34A":"#8c95c5", 
        "IDS1_KJEGr355l":"#8c5fa7", 
        "IDS1_rad":"#f26323"
    }

    compiled_df, filt_dfs = pa.master_histogram_func('test', output_folder=output_folder, thresh=FRET_thresh, swarmplot=False)
    compiled_df.groupby('treatment_name')['unique_id'].nunique()

    # --------------------------------- synchronised transitions ------------------------------------------
    
    compiled_df = pd.read_csv(f'{output_folder}/Cleaned_FRET_histogram_data.csv')
    order = list(compiled_df['treatment_name'].unique())
    list_to_drop = ['']
    filtered_list = [i for i in order if i not in list_to_drop]

    percent_trans_meet_criteria_df, calculated_transitions_df, consecutive_from_dnak_release, nonconsecutive_from_dnak_release, filt_data, sycnchronised_data, combined_consec_nonconsec = pg.master_plot_synchronised_transitions(order=order, 
                                        output_folder=output_folder, 
                                        exposure=0.2, 
                                        frames_to_plot=50, 
                                        FRET_before=FRET_thresh, 
                                        FRET_after=FRET_thresh, 
                                        datatype="Proportion", 
                                        filt=True, 
                                        palette=colors, 
                                        add_time=0) #### datatype could be ratio

    percent_trans_meet_criteria_df_sorted = percent_trans_meet_criteria_df.sort_values('treatment')
    percent_trans_meet_criteria_df_sorted.to_csv(f'{output_folder}/controlled_transitions.csv')

    ks_stats = pg.fit_ksdistribution(data=combined_consec_nonconsec, 
                                     output_folder=output_folder, 
                                     palette=colors)  #### palette is a dict of colors for each treatment

    fit_results = pg.extract_exp_kinetics(df=sycnchronised_data, 
                                       output_folder=output_folder, 
                                       exp_func=pg.exp_func, 
                                       x='time_from_trans', 
                                       y='FRET', 
                                       xlim=(0.5, 10), 
                                       palette=colors)  #### palette is a dict of colors for each treatment

    pg.plot_fit_parameters(fit_results, 
                           output_folder=output_folder,
                           palette=colors, 
                           list_to_keep=['IDS1_KJE', 'IDS1_KJEG', 'IDS1_rad'])

    # -------------------- combine timelapse data from different treatments and scripts -----------------------------------------------

    timepoint = {'6min':6,
                 '12min':12, 
                 '18min':18, 
                 '24min':24, 
                 '30min':30, 
                 '36min':36, 
                 '42min':42, 
                 '48min':48, 
                 '54min':54, 
                 '60min':60, 
                 }
    
    palette = {'IDS1':'BuPu', 
               'IDS2':'Blues', 
               'CDS1':'Oranges', 
               'NDS1':'Greens'}

    compiled_df_mapped, filt_dfs_mapped = pb.timelapse_mapping(output_folder, timepoint)
    pb.heatmap_timelapse(output_folder, compiled_df_mapped, palette)

    timelapse_colors = {
    'IDS1':'purple',    
    'IDS2':'skyblue', 
    'NDS1':'olivedrab',
    'CDS1':'darkorange',
    }

    # treatment_name in filt_dfs must be organized as such: 'protein_variable_timepoint'; where variable
    # is the experimental treatment (e.g., RNA variant or chaperone combination)
    filt_dfs, timepoint_plotdata, fit_dict = pb.master_timelapse_func(filt_dfs=filt_dfs_mapped, 
                             thresh=FRET_thresh, 
                             xlim_min=6, 
                             xlim_max=60, 
                             output_folder=output_folder, 
                             palette=timelapse_colors, 
                             data_type='normalised', 
                             markersize=75)  #### can either be raw values 'Mean' or normalised to 100% 'normalised'.

    rate_constants = pb.plot_rate_constants(output_folder, timelapse_colors, fit_dict)

    # ------------------- plot all traces with defined length or specific molecules of interest -------------------------------------
    
    data_paths_traces = {
        "example1":"trace1.dat",
        "example2":"trace2.dat",
    }

    pc.master_plot_individual_trace(data_paths_traces, input_folder=output_folder, exposure=0.2)
    pc.master_plot_traces_func(input_folder=output_folder, exposure=0.2, min_trace_length=200)

    # ------------------ process and plot data for live flow-in experiments -------------------------
    
    data_paths_liveflow = {
        'treatment':('treatment_name', 'Experiment_2-description/raw_data/treatment_to_plot'),
        'treatment':('treatment_name', 'Experiment_2-description/raw_data/treatment_to_plot'),
    }

    transition_type = 'low_to_high' ### or high_to_low

    pe.master_heatmap_processing_func(data_paths_liveflow, 
                                input_folder=output_folder, 
                                exposure=0.200, 
                                FRET_thresh=FRET_thresh, 
                                transition_type=transition_type, 
                                time_thresh=15, 
                                injection_time=15)

    pf.master_plot_flowin_func(input_folder=output_folder, 
                            transition_type=transition_type, 
                            gridsize=100, 
                            binshex=80)


    # --------------------------- gaussian fitting and plottting proportion of populations over time -------------------------------
    
    timepoint = {'6min':0, '12min':6,'18min':12,'24min':18,'30min':24,'36min':30}

    collated_df, gauss_saveloc = ph.fit_gauss_master(input_folder=output_folder, gauss_num='three', timepoint=None)
    test, collated_df, gauss_saveloc = ph.fit_gauss_master(input_folder=output_folder, gauss_num='three', timepoint=timepoint)

    # Combined the gaussian fit proportion data from multiple experiments and combine here to plot

    data_from_exp = {
    'KJE':'Experiment_1-KJE-timelapse/python_results/GaussianFits/collated_populations.csv',
    }

    final = ph.plot_gauss_timelapse(data_from_exp, gauss_saveloc)

    # --------------------------- FRET efficiency vs distance plot ----------------------------------

    r_values, fret_at_distance = pa.plot_fret_vs_distance(output_folder, 
                                                    R0=51, 
                                                    distance_range=(0, 120), 
                                                    num_points=500, 
                                                    fret_values=(0.05, 0.075), 
                                                    distance_to_mark=38)


#    # --------------------------- plot FRET heatmaps  -----------------------------------------------

    timepoint = {'0min':0,
                 '6min':6,
                 '12min':12, 
                 '18min':18, 
                 '24min':24, 
                 '30min':30, 
                 '36min':36, 
                 '42min':42, 
                 '48min':48, 
                 '54min':54, 
                 '60min':60, 
                 }
    
    protein_palette = {'IDS1':'BuPu', 
               'IDS2':'Blues', 
               'CDS1':'Oranges', 
               'NDS1':'Greens'}

    compiled_df_mapped, filt_dfs_mapped = pa.timelapse_mapping(output_folder, timepoint)

    pa.heatmap_timelapse(output_folder=output_folder, compiled_df=compiled_df_mapped, palette=protein_palette)
    compiled_df_mapped[compiled_df_mapped['treatment_name'].str.contains('native', case=False, na=False)]