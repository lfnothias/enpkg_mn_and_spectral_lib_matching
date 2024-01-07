import pandas as pd
import glob
import os
import math
import yaml
import shutil
import time
import sys
import logging
import textwrap
sys.path.append('src/')
import argparse

from matchms import Pipeline
from matchms.importing import load_from_mgf
from matchms.filtering import add_precursor_mz
from matchms.filtering.require_minimum_number_of_peaks  import require_minimum_number_of_peaks 

from spectral_db_loader import load_spectral_db
from spectral_db_loader import load_clean_spectral_db
from spectral_db_loader import save_spectral_db
from spectral_lib_matcher import spectral_matching
from molecular_networking import generate_mn
from install_ionmode_folder import check_for_ionmode_folder_and_restruct_if_needed

import math
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial
sys.stdout.flush()

print("Molecular networking and spectral library started", flush=True)

def load_configuration(config_path):
    """
    Load the configuration from a YAML file.
    
    Args:
        config_path (str): Path to the YAML configuration file.
    
    Returns:
        dict: Dictionary containing configuration parameters.
    """
    with open(config_path, 'r') as file:
        try:
            config = yaml.load(file, Loader=yaml.FullLoader)
        except yaml.YAMLError as e:
            print(f"Error loading configuration file: {e}")
            sys.exit(1)
        return config


def filter_spectra_by_ion_mode(spectra, ion_mode_start):
    filtered_spectra = []
    for spectrum in spectra:
        metadata = spectrum.metadata
        if metadata and 'ionmode' in metadata and metadata.get('ionmode', '').lower().startswith(ion_mode_start.lower()):
            filtered_spectra.append(spectrum)

    return filtered_spectra


# New function to load, filter, and save the spectral database
def load_filter_and_save_spectral_db(spectral_db_path, polarity, spectral_db_path_cleaned_polarity):
    print('Loading the spectral library: ' + spectral_db_path)
    spectral_db = load_clean_spectral_db(spectral_db_path)
    print('Number of spectra in the spectral library: ' + str(len(spectral_db)))

    print("Filtering by polarity")
    try:
        spectral_db = filter_spectra_by_ion_mode(spectral_db, polarity)
        print('Number of spectra in the spectral library in ' + polarity + ' mode: ' + str(len(spectral_db)))
        print('Saving the filtered spectral library: ' + spectral_db_path_cleaned_polarity)
        save_spectral_db(spectral_db, spectral_db_path_cleaned_polarity)
    except Exception as e:
        print(f"An error occurred: {e}")
        raise

    return spectral_db


def prepare_spectral_library(spectral_db_path, polarity):
    spectral_db_path_cleaned = spectral_db_path[:-4] + '_cleaned.msp'
    spectral_db_path_cleaned_polarity = spectral_db_path[:-4] + '_' + str(polarity) + '_cleaned.msp'

    if os.path.exists(spectral_db_path_cleaned_polarity):
        print('Already existing in polarity')
        print('Found a cleaned spectral library for the polarity already present: ' + spectral_db_path_cleaned_polarity)
        return load_clean_spectral_db(spectral_db_path_cleaned_polarity)
    elif os.path.exists(spectral_db_path_cleaned):
        print('Found a cleaned spectral library already present: ' + spectral_db_path_cleaned)
        return load_filter_and_save_spectral_db(spectral_db_path_cleaned, polarity, spectral_db_path_cleaned_polarity)
    else:
        print('Not cleaned yet')
        try:
            spectral_db_raw = load_spectral_db(spectral_db_path)
            print('Saving the cleaned spectral library: ' + spectral_db_path_cleaned)
            save_spectral_db(spectral_db_raw, spectral_db_path_cleaned)
            return load_filter_and_save_spectral_db(spectral_db_path_cleaned, polarity, spectral_db_path_cleaned_polarity)
        except Exception as e:
            print(f"An error occurred during spectral library preparation: {e}")
            raise


def copy_file_to_target(source_path, target_dir):
    """ Copy a file to a target directory. """
    try:
        if os.path.exists(source_path):
            target_path = os.path.join(target_dir, os.path.basename(source_path))
            shutil.copy(source_path, target_path)
            return target_path
        else:
            print(f"File not found: {source_path}")
            return None
    except Exception as e:
        print(f"Error copying file: {e}")
        return None

def process_sample(sample_dir, data_directory, polarity, mn_msms_mz_tol, mn_score_cutoff, mn_top_n, mn_max_links, recompute, spectral_db, parent_mz_tol, msms_mz_tol, min_score, min_peaks):
    try:
        script_dir = os.path.dirname(os.path.realpath(__file__))
        base_dir = os.path.dirname(script_dir)

        sample_dir_name = os.path.basename(sample_dir)

        config_file_path = os.path.join(base_dir, 'configs/user/user.yaml')
        params_list = load_configuration(config_file_path)

        if not recompute and os.path.isfile(lib_results_path):
            print(f"Skipping {sample_dir} as results already exist.")
            return sample_dir
    
        if not sample_dir_name:
            print("Sample directory name could not be extracted.")
            return

        spectra_file_path = glob.glob(os.path.join(data_directory, sample_dir, polarity, '*-FBMN.mgf'))
        feature_table_path = glob.glob(os.path.join(data_directory, sample_dir, polarity, '*-feature_table.csv'))

        # Correct path for spectra and feature table

        # Check for files in the parent directory and update paths if necessary
        if not spectra_file_path:
            spectra_file_path = glob.glob(os.path.join(data_directory, sample_dir, '*-FBMN.mgf'))
        if not feature_table_path:
            feature_table_path = glob.glob(os.path.join(data_directory, sample_dir, '*-feature_table.csv'))

        if feature_table_path:
            feature_table = pd.read_csv(feature_table_path[0], sep=',')
            # Additional processing code...
        else:
            print(f"No feature table found for {sample_dir}.")
            return

        print('''
        Treating files in the folder: ''' + sample_dir
        )

        mn_dir = os.path.join(data_directory, sample_dir, polarity, 'molecular_network')

        if not os.path.exists(mn_dir):
                os.makedirs(mn_dir)

        mn_ci_ouput_path = f'{mn_dir}/{sample_dir_name}_mn_metadata_'+polarity+'.tsv'
        mn_graphml_ouput_path = f'{mn_dir}/{sample_dir_name}_mn_'+polarity+'.graphml'
        mn_config_path = f'{mn_dir}/config.yaml'
        
        # Import query spectra
        spectra_query = list(load_from_mgf(spectra_file_path[0]))
        spectra_query = [require_minimum_number_of_peaks(s, n_required=1) for s in spectra_query]
        spectra_query = [add_precursor_mz(s) for s in spectra_query if s]

        # Molecular networking
        print('''
        ### Molecular networking 
        ''')
        start_time = time.time()

        if recompute:
        # If recompute is True, remove the content of mn_dir
            if os.path.exists(mn_dir):
                shutil.rmtree(mn_dir)
            os.makedirs(mn_dir)

        config_file_source_path = os.path.join(base_dir, 'configs/user/user.yaml')

        if not os.path.isfile(mn_graphml_ouput_path) or recompute:
            generate_mn(spectra_query, mn_graphml_ouput_path, mn_ci_ouput_path, mn_msms_mz_tol, mn_score_cutoff, mn_top_n, mn_max_links)
            shutil.copyfile(config_file_source_path, mn_config_path)

        end_time = time.time()
        elapsed_time = end_time - start_time
        print(f'''
        Molecular Networking done in {elapsed_time:.2f} secs.
        ''')

        print('''
        ### Spectral matching
        ''')
        start_time = time.time()

        lib_dir = os.path.join(data_directory, sample_dir, polarity, 'spectral_lib_matching')
        if not os.path.exists(lib_dir):
            os.makedirs(lib_dir)
        lib_results_path = f'{lib_dir}/{sample_dir_name}_lib_results_{polarity}.tsv'

        if recompute:
            # If recompute is True, remove the content of mn_dir
            print(f"Recomputing. Removing directory: {lib_dir}")
            shutil.rmtree(lib_dir)
            os.makedirs(lib_dir)

        print(f"Calling spectral_matching for: {lib_results_path}")
        try:
            spectral_matching(spectra_query, spectral_db, parent_mz_tol, msms_mz_tol, min_score, min_peaks, lib_results_path)
        except Exception as e:
            print(f"An error occurred during spectral matching: {e}")
            lib_results_path = None  # Set lib_results_path to None to indicate failure

        # Adding a short delay and checking if the file exists
        #time.sleep(1)  # wait for 1 seconds
        #if os.path.isfile(lib_results_path):
        #    continue
            #print(f"Confirmed creation of spectral matching results file: {lib_results_path}")
            #print("Absolute path of results file:", os.path.abspath(lib_results_path))
            #print(f"File size: {os.path.getsize(lib_results_path)} bytes")
        #else:
        #    print(f"File not found after creation attempt: {lib_results_path}")

        # Check if the spectral matching results file exists
        if os.path.isfile(lib_results_path):
            print(f"Spectral matching results file created: {lib_results_path}")
        else:
            print(f"No results file found at: {lib_results_path}")

        end_time = time.time()
        elapsed_time = end_time - start_time

        print(f'''
        Spectral library matching done in {elapsed_time:.2f} secs.
        ''')
        
        lib_dir = os.path.join(data_directory, sample_dir, polarity, 'spectral_lib_matching')

        #print(f"Checking existence of file: {lib_results_path}")
        #if os.path.exists(lib_results_path):
        #    print("File exists.")
        #else:
        #    print("File does not exist.")

        print("Reading dt_db_results from:", lib_results_path)
        try:
            dt_db_results = pd.read_csv(lib_results_path, sep='\t',
                                        usecols=['msms_score', 'feature_id', 'reference_id', 'adduct',
                                                'matched_peaks', 'charge',
                                                'ionmode', 'instrument', 'instrument_type', 'comment',
                                                'inchikey', 'inchi', 'smiles', 'compound_name', 'Spectral_library', 'Spectral_library_ID'], on_bad_lines='skip', low_memory=True)
        except FileNotFoundError as e:
            print(f"FileNotFoundError: {e}")
            print(f"Attempting to read from path: {lib_results_path}")
            raise
        except Exception as e:
            print(f"An unexpected error occurred: {e}")
            raise
            #dt_db_results = None  # Set dt_db_results to None to indicate failure

        print("Reading clusterinfo_summary from:", mn_ci_ouput_path)
        try:
            clusterinfo_summary = pd.read_csv(mn_ci_ouput_path, sep='\t',
                                            on_bad_lines='skip', low_memory=True)
            clusterinfo_summary.rename(columns={'precursor_mz': 'mz'}, inplace=True)
            #print("clusterinfo_summary loaded successfully")
        except Exception as e:
            print(f'Error reading clusterinfo_summary: {e}')
            clusterinfo_summary = None

        if dt_db_results is None:
            print("dt_db_results is None")
            print(dt_db_results)
        elif dt_db_results.empty:
            print("dt_db_results is empty")
        else:
            try:
                print("Preparing annotation results:")
                cluster_summary_path = f'{lib_dir}/{sample_dir_name}_cluster_info_{polarity}.tsv'
                clusterinfo_summary.to_csv(cluster_summary_path, sep='\t', index=False)

                #print("Before merging:")
                #print("dt_db_results shape:", dt_db_results.shape)
                #print("clusterinfo_summary shape:", clusterinfo_summary.shape)
                #print("clusterinfo_summary head:", clusterinfo_summary.head(5))
                dt_db_results = dt_db_results.sort_values(by='msms_score', ascending=False)  # Sort by 'msms_score' in descending order
                dt_db_results['feature_id'] = dt_db_results['feature_id'].astype(int)
                clusterinfo_summary['feature_id'] = clusterinfo_summary['feature_id'].astype(int)
                dt_db_results.reset_index(drop=True, inplace=True)
                clusterinfo_summary.reset_index(drop=True, inplace=True)
                # Merge and create the full table with all annotations
                dt_db_results_final = pd.merge(clusterinfo_summary, dt_db_results, on='feature_id', how='inner')
                spectral_db_path = os.path.join(base_dir, params_list['paths']['spectral_db_path'])
                dt_db_results_final['lib_name']= str(spectral_db_path)

                cols = ['feature_id'] + [col for col in dt_db_results_final.columns if col != 'feature_id']
                dt_db_results_final = dt_db_results_final[cols]

                # Save the full table
                db_results_full_path = f'{lib_dir}/{sample_dir_name}_lib_results_final_{polarity}.tsv'
                dt_db_results_final.to_csv(db_results_full_path, sep='\t', index=False)
                #print('Full results written out to:', db_results_full_path)

                # Drop duplicates for the unique table
                dt_db_results_sorted_unique = dt_db_results.sort_values(by='msms_score', ascending=False).drop_duplicates(subset='feature_id', keep='first')
                dt_db_results_sorted_unique.reset_index(drop=True, inplace=True)
                dt_db_results_unique = pd.merge(clusterinfo_summary, dt_db_results_sorted_unique, on='feature_id', how='inner')
                dt_db_results_unique['lib_name']= str(spectral_db_path)
                cols = ['feature_id'] + [col for col in dt_db_results_unique.columns if col != 'feature_id']
                dt_db_results_unique = dt_db_results_unique[cols]

                # Save the unique table
                db_results_unique_path = f'{lib_dir}/{sample_dir_name}_lib_results_final_{polarity}_unique.tsv'
                dt_db_results_unique.to_csv(db_results_unique_path, sep='\t', index=False)
                #print('Unique results written out to:', db_results_unique_path)

                # Print summary stats
                print('Number of features: ' + str(len(clusterinfo_summary)))
                print('Number of MS2 annotations in ' + polarity + ' mode: ' + str(len(dt_db_results)))
                print('Number of unique annotated features: ' + str(len(dt_db_results_unique['feature_id'].unique())))


                print('''
                ### Finished molecular networking and spectral library
                ''')
            except Exception as e:
                print(f"Unexpected error: {e}")
                raise  # Or handle it appropriately

    except Exception as e:
        print(f"Error processing sample {sample_dir}: {e}")
        return None  # Indicate failure

def process_samples(samples_dir, data_directory, polarity, spectral_db, mn_msms_mz_tol, mn_score_cutoff, mn_top_n, mn_max_links, recompute, parent_mz_tol, msms_mz_tol, min_score, min_peaks, num_cpus_to_use):
    valid_sample_dirs = []
    remaining = len(samples_dir)

    with ProcessPoolExecutor(max_workers=num_cpus_to_use) as executor:
        futures = {}
        for sample_dir in samples_dir:
            future = executor.submit(process_sample, sample_dir, data_directory, polarity, mn_msms_mz_tol, mn_score_cutoff, mn_top_n, mn_max_links, recompute, spectral_db, parent_mz_tol, msms_mz_tol, min_score, min_peaks)
            futures[future] = sample_dir

        for future in as_completed(futures):
            remaining -= 1
            print(f'Remaining folders: {remaining} out of {len(samples_dir)}')
            try:
                result = future.result()
                if result:
                    valid_sample_dirs.append(result)
            except Exception as e:
                print(f"A task failed with error: {e}")

    return valid_sample_dirs

def main():
    """ Argument parser """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent('''\
            This script runs molecular networking and spectral library matching on FBMN output files
                --------------------------------
                Arguments:
                - Path to the directory where samples folders are located
                - Ionization mode to process: pos / neg
            '''))

    parser.add_argument('-p', '--sample_dir_path', required=True,
                        help='The path to the directory where samples folders to process are located')
    parser.add_argument('-ion', '--ionization_mode', required=True,
                        choices=['pos', 'neg'],
                        help='The ionization mode to perform spectral library matching')
    parser.add_argument('-c', '--cpus', type=int, default=None,
                    help='Number of CPUs to use for processing. Default is 80% of available CPUs.')
    parser.add_argument('-r', '--recompute', action='store_true', 
                        help='Recompute even if the files are already present')
    

    pd.options.mode.chained_assignment = None

    os.chdir(os.getcwd())

    args = parser.parse_args()
    data_directory = os.path.normpath(args.sample_dir_path)
    polarity = args.ionization_mode
    recompute = args.recompute

    if len(sys.argv) < 3:
        print("""
              Usage: 
              python enpkg_mn_and_matching.py -p <data_directory> -ion <polarity>
              python enpkg_mn_and_matching.py -p <data_directory> -ion <polarity>  -c <cpu to use> -r <recompute>
              """)
        sys.exit(1)

    ###### START #####
    print('=== Start of ENPKG Molecular networking and spectral library matching ====')
    print('    ')

    script_dir = os.path.dirname(os.path.realpath(__file__))
    base_dir = os.path.dirname(script_dir)
    config_file_path = os.path.join(base_dir, 'configs/user/user.yaml')
    params_list = load_configuration(config_file_path)

    spectral_db_path = os.path.join(base_dir, params_list['paths']['spectral_db_path'])
    print('> Preparing the spectral library:')
    print('  This is done only once entirely')
    print('Found a spectral library: '+spectral_db_path)
    spectral_db = prepare_spectral_library(spectral_db_path, polarity)

    # Sample processing
    check_for_ionmode_folder_and_restruct_if_needed(data_directory, polarity)


    # Parallel processing
    # Get user-specified number of CPUs or use default
    total_cpus = os.cpu_count()
    print('Total number of cpus: ',  total_cpus)
    user_cpus = args.cpus
    num_cpus_to_use = user_cpus if user_cpus else math.ceil(total_cpus * 0.8)

    print(f'Number of CPUs used: {num_cpus_to_use}')

    # Rest of the configuration parameters
    recompute = params_list['general_params']['recompute']

    samples_dir = [os.path.join(sample_dir) for sample_dir in os.listdir(data_directory) if os.path.isdir(os.path.join(data_directory, sample_dir))]

    print(f'{len(samples_dir)} samples folder were detected in the input directory. They will be checked for minimal requirements.')

    valid_sample_dirs = []
    
    # Minimal checks
    for sample_dir in samples_dir:
        full_sample_path = os.path.join(data_directory, sample_dir)
        polarity_path = os.path.join(full_sample_path, polarity)

        # Ensure the polarity subfolder exists
        if not os.path.exists(polarity_path):
            os.makedirs(polarity_path)

        # Check for files in the polarity subfolder
        spectra_file_paths = glob.glob(os.path.join(polarity_path, '*-FBMN.mgf'))
        feature_table_paths = glob.glob(os.path.join(polarity_path, '*-feature_table.csv'))

        # If files not found in the polarity subfolder, check the parent directory and copy
        if not spectra_file_paths:
            spectra_file_paths_up = glob.glob(os.path.join(full_sample_path, '*-FBMN.mgf'))
            for file_path in spectra_file_paths_up:
                copied_path = copy_file_to_target(file_path, polarity_path)
                if copied_path:
                    spectra_file_paths.append(copied_path)

        if not feature_table_paths:
            feature_table_paths_up = glob.glob(os.path.join(full_sample_path, '*-feature_table.csv'))
            for file_path in feature_table_paths_up:
                copied_path = copy_file_to_target(file_path, polarity_path)
                if copied_path:
                    feature_table_paths.append(copied_path)

        # Check again after copying
        if not spectra_file_paths:
            print(f"{full_sample_path} folder has no MSMS data, it is removed from the processing list.")
            continue

        if not feature_table_paths:
            print(f"{full_sample_path} folder has no feature intensity table, it is removed from the processing list.")
            continue

        valid_sample_dirs.append(full_sample_path)

    print(f'{len(valid_sample_dirs)} samples folder were found to be complete and will be processed.')

    process_samples(valid_sample_dirs, data_directory, polarity, spectral_db, 
                                        params_list['networking_params']['mn_msms_mz_tol'], 
                                        params_list['networking_params']['mn_score_cutoff'], 
                                        params_list['networking_params']['mn_top_n'], 
                                        params_list['networking_params']['mn_max_links'], 
                                        recompute, 
                                        params_list['spectral_match_params']['parent_mz_tol'], 
                                        params_list['spectral_match_params']['msms_mz_tol'], 
                                        params_list['spectral_match_params']['min_score'], 
                                        params_list['spectral_match_params']['min_peaks'],
                                        num_cpus_to_use)
    
    
    print("All tasks have been processed successfully for the following!:", valid_sample_dirs)

if __name__ == '__main__':
    main()