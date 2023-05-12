import pandas as pd
import glob
import os
import yaml
import shutil
import time
import sys
import logging
sys.path.append('src/')

from matchms import Pipeline
from matchms.importing import load_from_mgf
from matchms.filtering import add_precursor_mz
from matchms.filtering.require_minimum_number_of_peaks  import require_minimum_number_of_peaks 

from spectral_db_loader import load_spectral_db
from spectral_db_loader import load_clean_spectral_db
from spectral_db_loader import save_spectral_db
from spectral_lib_matcher import spectral_matching
from molecular_networking import generate_mn

pd.options.mode.chained_assignment = None

os.chdir(os.getcwd())

if len(sys.argv) < 2:
    print("Usage: python enpkg_mn_and_matching.py <data_directory> ")
    sys.exit(1)

data_directory = sys.argv[1]

# Rest of the code remains unchanged...


# you can copy the configs/default/default.yaml to configs/user/user.yaml

with open (r'configs/user/user.yaml') as file:    
    params_list = yaml.load(file, Loader=yaml.FullLoader)

recompute = params_list['general_params']['recompute']

spectral_db_path = params_list['paths']['spectral_db_path']

parent_mz_tol = params_list['spectral_match_params']['parent_mz_tol']
msms_mz_tol = params_list['spectral_match_params']['msms_mz_tol']
min_score = params_list['spectral_match_params']['min_score']
min_peaks = params_list['spectral_match_params']['min_peaks']

mn_msms_mz_tol = params_list['networking_params']['mn_msms_mz_tol']
mn_score_cutoff = params_list['networking_params']['mn_score_cutoff']
mn_max_links = params_list['networking_params']['mn_max_links']
mn_top_n = params_list['networking_params']['mn_top_n']

###### START #####
print('=== Start of ENPKG Molecular networking and spectral library matching ====')
# Prepare spectral library (only once)
print('    ')
print('> Preparing the spectral library')
print('  This is done only once')
print('Found a spectral library: '+spectral_db_path)
spectral_db_path_cleaned = spectral_db_path[:-4] + '_cleaned.msp'
if os.path.exists(spectral_db_path_cleaned):
    start_time = time.time()
    print('Found a cleaned spectral library already present: '+spectral_db_path_cleaned)
    ""
    spectral_db = load_clean_spectral_db(spectral_db_path_cleaned)

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f'''
    Loading in {elapsed_time:.2f} secs.
    ''')
else:
    start_time = time.time()
    print('Loading the spectral library: '+spectral_db_path)
    spectral_db_raw = load_spectral_db(spectral_db_path)
    print('Saving the cleaned spectral library: '+spectral_db_path_cleaned)
    save_spectral_db(spectral_db_raw, spectral_db_path_cleaned)
    print('Loading the cleaned spectral library: '+spectral_db_path)
    spectral_db = load_clean_spectral_db(spectral_db_path_cleaned)
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f'''
    Cleaning and loading in {elapsed_time:.2f} secs.
    ''')


# Iteration over samples directory to count samples with required input files

samples_dir = [directory for directory in os.listdir(data_directory)]
print(f'{len(samples_dir)} samples folder were detected in the input directory. They will be checked for minimal requirements.')

for sample_dir in samples_dir:
    print(sample_dir)
    # Construct path using os.path.join()
    spectra_file_path = glob.glob(os.path.join(data_directory, sample_dir, '*.mgf'))

    feature_table_path = glob.glob(os.path.join(data_directory, sample_dir, '*feature_table.csv'))
    # Check if MS/MS spectra are present 
    if len(spectra_file_path) != 0 :
        pass
    else:
        print(sample_dir + "folder has no MSMS data, it is removed from the processing list.")
        samples_dir.remove(sample_dir)
        continue

    if len(feature_table_path) != 0 :
        pass
    else:
        print(sample_dir + "folder  has no feature intensity table, it is removed from the processing list.")
        samples_dir.remove(sample_dir)
        continue

    # else:
    #     continue

print(f'{len(samples_dir)} samples folder were found to be complete and will be processed.')

for sample_dir in samples_dir:
    print(sample_dir)
    results_path = os.path.join(data_directory, sample_dir)
    if not os.path.exists(results_path):
        os.makedirs(results_path)

    spectra_file_path = glob.glob(os.path.join(data_directory, sample_dir, '*.mgf'))[0]
    feature_table_path = glob.glob(os.path.join(data_directory, sample_dir, '*feature_table.csv'))[0]
    print(feature_table_path)
    feature_table = pd.read_csv(feature_table_path, sep=',')
    # Additional processing code goes here

    print('''
    Treating files in the folder: ''' + sample_dir
    )

    mn_ci_ouput_path = f'{results_path}/molecular_network/{sample_dir}_mn_metadata.tsv'
    mn_graphml_ouput_path = f'{results_path}/molecular_network/{sample_dir}_mn.graphml'
    mn_config_path = f'{results_path}/molecular_network/config.yaml'
    mn_folder_path = f'{results_path}/molecular_network/'
    
    # Import query spectra
    spectra_query = list(load_from_mgf(spectra_file_path))
    spectra_query = [require_minimum_number_of_peaks(s, n_required=1) for s in spectra_query]
    spectra_query = [add_precursor_mz(s) for s in spectra_query if s]

    # Molecular networking
    print('''
    Molecular networking 
    ''')
    start_time = time.time()

    generate_mn(spectra_query, mn_graphml_ouput_path, mn_ci_ouput_path, mn_msms_mz_tol, mn_score_cutoff, mn_top_n, mn_max_links)
    shutil.copyfile(r'configs/user/user.yaml', mn_config_path)
    
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f'''
    Molecular Networking done in {elapsed_time:.2f} secs.
    ''')

    print('''
    Spectral matching
    ''')
    start_time = time.time()
    
    lib_results_path  = f'{results_path}/spectral_lib_matching/{sample_dir}_lib_results.tsv'

    print(spectra_query)
    spectral_matching(spectra_query, spectral_db, parent_mz_tol,
        msms_mz_tol, min_score, min_peaks, lib_results_path)
    
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f'''
    Spectral matching done in {elapsed_time:.2f} secs.
    ''')


    try:
        dt_db_results = pd.read_csv(lib_results_path, sep='\t',
            usecols=['msms_score', 'feature_id', 'reference_id', 'short_inchikey'], on_bad_lines='skip', low_memory=True)
    except:   
        continue
    # Add 'libname' column and rename msms_score column
    dt_db_results['libname'] = 'PublicSpecLibrary'
    # Load MN metadata
    clusterinfo_summary = pd.read_csv(mn_ci_ouput_path, sep='\t', usecols=['feature_id', 'precursor_mz', 'component_id'], \
        on_bad_lines='skip', low_memory=True)
    clusterinfo_summary.rename(columns={'precursor_mz': 'mz'}, inplace=True)

    clusterinfo_summary.to_csv(lib_results_path[:-4]+'_cluster_info.tsv', sep='\t', index=False)
    dt_db_results_final = pd.merge(dt_db_results, clusterinfo_summary, on='feature_id')

    # Reorder the columns with 'feature_id' as the first column
    cols = ['feature_id'] + [col for col in dt_db_results_final.columns if col != 'feature_id']
    dt_db_results_final = dt_db_results_final[cols]

    dt_db_results_final.to_csv(lib_results_path[:-4]+'_final.tsv', sep='\t', index=False)


    print('Number of features: ' + str(len(clusterinfo_summary)))
    print('Number of MS2 annotation: ' + str(len(dt_db_results)))
    print('Number of annotated features: ' + str(len(dt_db_results['feature_id'].unique())))

             
    print('''
    Finished file: ''' + sample_dir
    )