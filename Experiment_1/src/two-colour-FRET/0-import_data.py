from smfret.src.processing_scripts import data_import_script as ps

rawdata_folder = "Experiment_1-description/raw_data"  ## Change for each experiment
input_folders = ps.locate_raw_drive_files(f'{rawdata_folder}/raw_data.txt')
ps.move_folders(input_folders, '.dat', rawdata_folder)

