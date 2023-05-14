import pandas as pd
import glob
import os
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

""" Argument parser """
parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description=textwrap.dedent('''\
        This script generate a RDF graph (.ttl format) from samples' individual Sirius annotations
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

pd.options.mode.chained_assignment = None

os.chdir(os.getcwd())

args = parser.parse_args()
data_directory = os.path.normpath(args.sample_dir_path)
polarity = args.ionization_mode

if len(sys.argv) < 3:
    print("Usage: python enpkg_mn_and_matching.py -p <data_directory> -ion <polarity>")
    sys.exit(1)

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
print('  This is done only once entirely')
print('Found a spectral library: '+spectral_db_path)
spectral_db_path_cleaned = spectral_db_path[:-4] +'_cleaned.msp'
spectral_db_path_cleaned_polarity = spectral_db_path[:-4] + '_'+str(polarity)+'_cleaned.msp'

def filter_spectra_by_ion_mode(spectra, ion_mode_start):
    filtered_spectra = []

    for spectrum in spectra:
        if spectrum.metadata.get('ionmode', '').lower().startswith(ion_mode_start.lower()):
            filtered_spectra.append(spectrum)

    return filtered_spectra

if os.path.exists(spectral_db_path_cleaned_polarity):
    print('Already existing in polarity')
    print('Found a cleaned spectral library for the polarity already present: '+spectral_db_path_cleaned_polarity)
    spectral_db = load_clean_spectral_db(spectral_db_path_cleaned_polarity)
    print(' Number of spectra in the spectral library in '+polarity+' mode:' +str(len(spectral_db)))

elif os.path.exists(spectral_db_path_cleaned):
    start_time = time.time()
    print('Found a cleaned spectral library already present: '+spectral_db_path_cleaned)

    spectral_db = load_clean_spectral_db(spectral_db_path_cleaned)
    print(' Number of spectra in the spectral library in '+str(len(spectral_db)))
    
    print("Filtering by polarity")
    try:
        spectral_db = filter_spectra_by_ion_mode(spectral_db, polarity)
        print(' Number of spectra in the spectral library in '+polarity+' mode:' +str(len(spectral_db)))
        print('Saving the filtered spectral library: '+spectral_db_path_cleaned_polarity)
        save_spectral_db(spectral_db, spectral_db_path_cleaned_polarity)

    except Exception as e:
        print(f"An error occurred: {e}")

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f'''
    Loading in {elapsed_time:.2f} secs.
    ''')
else:
    print('Not cleaned yet')
    start_time = time.time()
    print('Loading the spectral library: '+spectral_db_path)
    spectral_db_raw = load_spectral_db(spectral_db_path)

    print('Saving the cleaned spectral library: '+spectral_db_path_cleaned)
    save_spectral_db(spectral_db_raw, spectral_db_path_cleaned)

    print('Loading the cleaned spectral library: '+spectral_db_path)
    spectral_db = load_clean_spectral_db(spectral_db_path_cleaned)

    print(' Number of spectra in the spectral library in '+str(len(spectral_db)))
    
    print("Filtering by polarity")
    try:
        spectral_db = filter_spectra_by_ion_mode(spectral_db, polarity)
        print(' Number of spectra in the spectral library in '+polarity+' mode:' +str(len(spectral_db)))
        print('Saving the filtered spectral library: '+spectral_db_path_cleaned_polarity)
        save_spectral_db(spectral_db, spectral_db_path_cleaned_polarity)
    except Exception as e:
        print(f"An error occurred: {e}")

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f'''
    Loading in {elapsed_time:.2f} secs.
    ''')

# Iteration over samples directory to count samples with required input files

check_for_ionmode_folder_and_restruct_if_needed(data_directory, polarity)

samples_dir = [os.path.join(sample_dir, polarity) 
               for sample_dir in os.listdir(data_directory)
               if os.path.isdir(os.path.join(data_directory, sample_dir))]

print(samples_dir)

print(f'{len(samples_dir)} samples folder were detected in the input directory. They will be checked for minimal requirements.')

valid_sample_dirs = []

for sample_dir in samples_dir:
    spectra_file_path = glob.glob(os.path.join(data_directory, sample_dir, '*ms2_'+polarity+'.mgf'))
    feature_table_path = glob.glob(os.path.join(data_directory, sample_dir, '*quant_'+polarity+'.csv'))

    if len(spectra_file_path) == 0 :
        print(sample_dir + " folder has no MSMS data, it is removed from the processing list.")
        continue

    if len(feature_table_path) == 0 :
        print(sample_dir + " folder  has no feature intensity table, it is removed from the processing list.")
        continue

    valid_sample_dirs.append(sample_dir)

samples_dir = valid_sample_dirs

print(f'{len(samples_dir)} samples folder were found to be complete and will be processed.')

for sample_dir in samples_dir:

    sample_dir_name = sample_dir.split('/')
    if len(sample_dir_name) > 1:
        sample_dir_name = sample_dir_name[0]
    else:
        print("No '/' character found in the string.")

    spectra_file_path = glob.glob(os.path.join(data_directory, sample_dir, '*ms2_'+polarity+'.mgf'))[0]
    feature_table_path = glob.glob(os.path.join(data_directory, sample_dir, '*quant_'+polarity+'.csv'))[0]
    feature_table = pd.read_csv(feature_table_path, sep=',')
    # Additional processing code goes here

    print('''
    Treating files in the folder: ''' + sample_dir
    )

    mn_dir = os.path.join(data_directory, sample_dir, 'molecular_network')
    if not os.path.exists(mn_dir):
            os.makedirs(mn_dir)

    mn_ci_ouput_path = f'{mn_dir}/{sample_dir_name}_mn_metadata_'+polarity+'.tsv'
    mn_graphml_ouput_path = f'{mn_dir}/{sample_dir_name}_mn_'+polarity+'.graphml'
    mn_config_path = f'{mn_dir}/config.yaml'
    
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
    
    lib_dir = os.path.join(data_directory, sample_dir, 'spectral_lib_matching')
    if not os.path.exists(mn_dir):
            os.makedirs(mn_dir)

    lib_results_path  = f'{lib_dir}/{sample_dir_name}_lib_results_'+polarity+'.tsv'

    spectral_matching(spectra_query, spectral_db, parent_mz_tol,
        msms_mz_tol, min_score, min_peaks, lib_results_path)
    
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f'''
    Spectral matching done in {elapsed_time:.2f} secs.
    ''')

    try:
        dt_db_results = pd.read_csv(lib_results_path, sep='\t',
            usecols=['msms_score', 'feature_id', 'reference_id',
                     'matched_peaks', 'adduct',	'charge',
                     'ionmode',	'instrument',	'instrument_type',	'comment',
                     'inchikey', 'inchi',	'smiles','compound_name'], on_bad_lines='skip', low_memory=True)
    except:   
        continue

    # Add 'libname' column and rename msms_score column
    dt_db_results['libname'] = spectral_db_path
    # Load MN metadata
    clusterinfo_summary = pd.read_csv(mn_ci_ouput_path, sep='\t', \
        on_bad_lines='skip', low_memory=True)
    clusterinfo_summary.rename(columns={'precursor_mz': 'mz'}, inplace=True)

    clusterinfo_summary.to_csv(lib_results_path[:-4]+'_cluster_info_'+polarity+'.tsv', sep='\t', index=False)
    dt_db_results_final = pd.merge(dt_db_results, clusterinfo_summary, on='feature_id')

    # Reorder the columns with 'feature_id' as the first column
    cols = ['feature_id'] + [col for col in dt_db_results_final.columns if col != 'feature_id']
    dt_db_results_final = dt_db_results_final[cols]

    dt_db_results_final.to_csv(lib_results_path[:-4]+'_final_'+polarity+'.tsv', sep='\t', index=False)

    print('Number of features: ' + str(len(clusterinfo_summary)))
    print('Number of MS2 annotation in '+polarity+' mode: ' + str(len(dt_db_results)))
    print('Number of annotated features: ' + str(len(dt_db_results['feature_id'].unique())))

             
    print('''
    Finished molecular networking and spectral library
    ''')