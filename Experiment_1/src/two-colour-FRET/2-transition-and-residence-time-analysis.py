from smfret.src.plotting_scripts import TDP_plots as pa
from smfret.src.plotting_scripts import dwelltime_plots as pb
from smfret.src.plotting_scripts import dwelltime_plots as pc
from smfret.src.plotting_scripts import residence_time_plots as pe

if __name__ == "__main__":

    output_folder = 'Experiment_1-description/python_results'
    order = ['treatment', 'treatment2', 'treatment3']

    # ------------------------------- generate transitions and plot TDPs ----------------------------------------------------------------
    pa.master_tdp_cleanup_func(output_folder=output_folder, 
                            exposure=0.2, 
                            FRET_value_for_classes=0.3,
                            FRET_to_filt_tdp=0.3)


    pa.master_TDP_plot(input_folder=output_folder, filt=True)


    # ----------------------------- plot dwell time data (e.g. bind-and-release, FRET before or after transition) -----------------------


    bind_release_col = pb.master_dwell_time_func(order=order, 
                        output_folder=output_folder, 
                        palette_to_use='BuPu', 
                        FRET_thresh=0.3, #### FRET value at which to filter data above or below. 
                        thresh_for_events=0.3, #### FRET value to threshold for count chaperone events (i.e., binding and unbinding)
                        fps=5,  ### frames per second
                        thresh=2, ### should be 10x expsoure if using NL515 smoothing on MASH FRET
                        Transition_threshold=0.3,
                        event='binding_and_release'
    )     


    # ----------------------------- plot transition frequence data ---------------------------------------------------------------

    sorted_list = sorted(order)   
    index_mapping = {value: index for index, value in enumerate(sorted_list)}
    indexes = [index_mapping[value] for value in order]
    pc.transition_frequency_plot(output_folder, indexes, FRET_thresh=0.5)

    # ---------------------------- plot residence time data (reports as mean +- SE, violin and cumulative histogram) -------------------------------------------
    
    data_paths_violin = {}
    for treatment in order:
        data_paths_violin[treatment] = f"{output_folder}/Dwell_times/Filtered_dwelltime_{treatment}.csv"

    colors_violin = {
        'RNA01': 'black', 
        'RNA09':"#ee9b00", 
        'RNA10':"#ca6702", 
        'RNA22':"#ae2012", 
        'RNA22duplex':"royalblue",
        # 'DNA02':"grey"
    }

    final, collated, cumulative_dwell_filt, halftime_summary = pe.master_residence_time_func(output_folder, 
                                                                        data_paths_violin, 
                                                                        order,
                                                                        palette='BuPu', 
                                                                        FRET_thresh=0.3, 
                                                                        binwidth=10, 
                                                                        cumulative_hist_binwidth=2, 
                                                                        fit_xlim=300, 
                                                                        plot_xlim=100)
