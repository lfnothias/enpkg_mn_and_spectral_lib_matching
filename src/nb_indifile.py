import pandas as pd
import glob
import os
import yaml
import shutil

from matchms.importing import load_from_mgf
from matchms.filtering import add_precursor_mz
from matchms.filtering.require_minimum_number_of_peaks  import require_minimum_number_of_peaks 

from spectral_db_loader import load_spectral_db
from spectral_lib_matcher import spectral_matching
from molecular_networking import generate_mn
from ms1_matcher import ms1_matcher
from reweighting_functions import taxonomical_reponderator, chemical_reponderator
from helpers import top_N_slicer, annotation_table_formatter

from plotter import plotter_count, plotter_intensity
from formatters import feature_intensity_table_formatter

pd.options.mode.chained_assignment = None

os.chdir(os.getcwd())

# for debug ony should be commented later
from pathlib import Path
p = Path(__file__).parents[1]
os.chdir(p)


# you can copy the configs/default/default.yaml to configs/user/user.yaml

with open (r'configs/user/user.yaml') as file:    
    params_list = yaml.load(file, Loader=yaml.FullLoader)

recompute = params_list['general_params'][0]['recompute']
ionization_mode = params_list['general_params'][1]['ionization_mode']

repository_path = params_list['paths'][0]['repository_path']
taxo_db_metadata_path = params_list['paths'][1]['taxo_db_metadata_path']
spectral_db_path = params_list['paths'][2]['spectral_db_path']
adducts_pos_path = params_list['paths'][3]['adducts_pos_path']
adducts_neg_path = params_list['paths'][4]['adducts_neg_path']

parent_mz_tol = params_list['spectral_match_params'][0]['parent_mz_tol']
msms_mz_tol = params_list['spectral_match_params'][1]['msms_mz_tol']
min_score = params_list['spectral_match_params'][2]['min_score']
min_peaks = params_list['spectral_match_params'][3]['min_peaks']

mn_parent_mz_tol = params_list['networking_params'][0]['mn_parent_mz_tol']
mn_msms_mz_tol = params_list['networking_params'][1]['mn_msms_mz_tol']
mn_score_cutoff = params_list['networking_params'][2]['mn_score_cutoff']
mn_max_links = params_list['networking_params'][3]['mn_max_links']
mn_top_n = params_list['networking_params'][4]['mn_top_n']

top_to_output= params_list['reweighting_params'][0]['top_to_output']
ppm_tol_ms1 = params_list['reweighting_params'][1]['ppm_tol_ms1']
use_post_taxo = params_list['reweighting_params'][2]['use_post_taxo']
top_N_chemical_consistency = params_list['reweighting_params'][3]['top_N_chemical_consistency']
min_score_taxo_ms1 = params_list['reweighting_params'][4]['min_score_taxo_ms1']
min_score_chemo_ms1 = params_list['reweighting_params'][5]['min_score_chemo_ms1']


###### START #####

# Iteration over samples directory to count samples with required input files

samples_dir = [directory for directory in os.listdir(repository_path)]

i = 0 
for sample_dir in samples_dir:
    try:
        metadata_file_path = os.path.join(repository_path, sample_dir, sample_dir + '_metadata.tsv')
    except FileNotFoundError:
        continue
    # Check if MS/MS spectra are present 
    if len(glob.glob(repository_path + sample_dir + '/' + ionization_mode + '/*'+ '_features_ms2_' + ionization_mode + '.mgf')) != 0 :
        spectra_file_path = glob.glob(repository_path + sample_dir + '/' + ionization_mode + '/*'+ '_features_ms2_' + ionization_mode + '.mgf')[0]
    
    # Check if features intensity table is present
    if len(glob.glob(repository_path + sample_dir + '/' + ionization_mode + '/*'+ '_features_quant_' + ionization_mode + '.csv')) != 0 :
        feature_table_path = glob.glob(repository_path + sample_dir + '/' + ionization_mode + '/*'+ '_features_quant_' + ionization_mode + '.csv')[0]
        feature_table = pd.read_csv(feature_table_path, sep=',')
    
    # Check if taxonomical metadata are present
    if len(glob.glob(repository_path + sample_dir + '/taxo_output' + '/*'+ '_taxo_metadata.tsv')) != 0 :
        taxo_metadata_path = glob.glob(repository_path + sample_dir + '/taxo_output' + '/*'+ '_taxo_metadata.tsv')[0]
        taxo_metadata = pd.read_csv(taxo_metadata_path, sep='\t')        
    else:
        continue
    i += 1
    
print(f'{i} samples with required input files detected')
if input("Do you wish to continue and process samples? (y/n)") != ("y"):
    exit()
    
# Load spectral DB
spectral_db = load_spectral_db(spectral_db_path)

# Calculate min and max m/z value using user's tolerance for adducts search
if ionization_mode == 'pos':
    adducts_df = pd.read_csv(adducts_pos_path, compression='gzip', sep='\t')
elif ionization_mode == 'neg':
    adducts_df = pd.read_csv(adducts_neg_path, compression='gzip', sep='\t')
else:
    raise ValueError('ionization_mode parameter must be pos or neg')

adducts_df['min'] = adducts_df['adduct_mass'] - \
    int(ppm_tol_ms1) * (adducts_df['adduct_mass'] / 1000000)
adducts_df['max'] = adducts_df['adduct_mass'] + \
    int(ppm_tol_ms1) * (adducts_df['adduct_mass'] / 1000000)

# Load structures taxonimcal data
if taxo_db_metadata_path.endswith('.csv.gz'):
    db_metadata = pd.read_csv(taxo_db_metadata_path, sep=',', compression='gzip', error_bad_lines=False, low_memory=False)
elif taxo_db_metadata_path.endswith('.csv'):
    db_metadata = pd.read_csv(taxo_db_metadata_path, sep=',', error_bad_lines=False, low_memory=False)
db_metadata['short_inchikey'] = db_metadata.structure_inchikey.str.split(
    "-", expand=True)[0]
db_metadata.reset_index(inplace=True)
    
# Processing
for sample_dir in samples_dir:
    try:
        metadata_file_path = os.path.join(repository_path, sample_dir, sample_dir + '_metadata.tsv')
    except FileNotFoundError:
        continue
     
    if len(glob.glob(repository_path + sample_dir + '/' + ionization_mode + '/*'+ '_features_ms2_' + ionization_mode + '.mgf')) != 0 :
        spectra_file_path = glob.glob(repository_path + sample_dir + '/' + ionization_mode + '/*'+ '_features_ms2_' + ionization_mode + '.mgf')[0]
        
    if len(glob.glob(repository_path + sample_dir + '/' + ionization_mode + '/*'+ '_features_quant_' + ionization_mode + '.csv')) != 0 :
        feature_table_path = glob.glob(repository_path + sample_dir + '/' + ionization_mode + '/*'+ '_features_quant_' + ionization_mode + '.csv')[0]
        feature_table = pd.read_csv(feature_table_path, sep=',')
        
    if len(glob.glob(repository_path + sample_dir + '/taxo_output' + '/*'+ '_taxo_metadata.tsv')) != 0 :
        taxo_metadata_path = glob.glob(repository_path + sample_dir + '/taxo_output' + '/*'+ '_taxo_metadata.tsv')[0]
        taxo_metadata = pd.read_csv(taxo_metadata_path, sep='\t')       
    else:
        continue

    print('''
    Treating file: ''' + sample_dir
    )

    isdb_results_path = f'{repository_path}{sample_dir}/{ionization_mode}/isdb/{sample_dir}_isdb_{ionization_mode}.tsv'
    mn_ci_ouput_path = f'{repository_path}{sample_dir}/{ionization_mode}/molecular_network/{sample_dir}_mn_metadata_{ionization_mode}.tsv'
    repond_table_path = f'{repository_path}{sample_dir}/{ionization_mode}/isdb/{sample_dir}_isdb_reweighted_{ionization_mode}.tsv'
    repond_table_flat_path = f'{repository_path}{sample_dir}/{ionization_mode}/isdb/{sample_dir}_isdb_reweighted_flat_{ionization_mode}.tsv'
    mn_graphml_ouput_path = f'{repository_path}{sample_dir}/{ionization_mode}/molecular_network/{sample_dir}_mn_{ionization_mode}.graphml'
    treemap_chemo_counted_results_path = f'{repository_path}{sample_dir}/{ionization_mode}/isdb/{sample_dir}_treemap_chemo_counted_{ionization_mode}.html'
    treemap_chemo_intensity_results_path = f'{repository_path}{sample_dir}/{ionization_mode}/isdb/{sample_dir}_treemap_chemo_intensity{ionization_mode}.html'
    isdb_config_path = f'{repository_path}{sample_dir}/{ionization_mode}/isdb/config.yaml'
    mn_config_path = f'{repository_path}{sample_dir}/{ionization_mode}/molecular_network/config.yaml'

    # Import query spectra
    spectra_query = list(load_from_mgf(spectra_file_path))
    spectra_query = [require_minimum_number_of_peaks(s, n_required=1) for s in spectra_query]
    spectra_query = [add_precursor_mz(s) for s in spectra_query if s]
    
    # Spectral matching
    print('''
    Spectral matching
    ''')
    
    spectral_matching(spectra_query, spectral_db, parent_mz_tol,
        msms_mz_tol, min_score, min_peaks, isdb_results_path)
    
    print('''
    Spectral matching done
    ''')
    
    # Molecular networking
    print('''
    Molecular networking 
    ''')
    
    generate_mn(spectra_query, mn_graphml_ouput_path, mn_ci_ouput_path, mn_msms_mz_tol, mn_score_cutoff, mn_top_n, mn_max_links)

    print('''
    Molecular Networking done
    ''')

    # Load ISDB results
    dt_isdb_results = pd.read_csv(isdb_results_path, sep='\t', \
        usecols=['msms_score', 'feature_id', 'reference_id', 'short_inchikey'], error_bad_lines=False, low_memory=True)
    
    # Add 'libname' column and rename msms_score column
    dt_isdb_results['libname'] = 'ISDB'

    # Load MN metadata
    clusterinfo_summary = pd.read_csv(mn_ci_ouput_path, sep='\t', usecols=['feature_id', 'precursor_mz', 'component_id'], \
        error_bad_lines=False, low_memory=True)
    clusterinfo_summary.rename(columns={'precursor_mz': 'mz'}, inplace=True)

    # Merge this back with the isdb results 
    dt_isdb_results = pd.merge(dt_isdb_results, clusterinfo_summary, on='feature_id')

    print('Number of features: ' + str(len(clusterinfo_summary)))
    print('Number of MS2 annotation: ' + str(len(dt_isdb_results)))
    print('Number of annotated features: ' + str(len(dt_isdb_results['feature_id'].unique())))
    
    # Now we directly do the MS1 matching stage on the cluster_summary. 
    print('''
    MS1 annotation
    ''')
    
    df_MS1 = ms1_matcher(clusterinfo_summary, adducts_df, db_metadata)
    
    print('''
    MS1 annotation done
    ''')
     
    # Merge MS1 results with MS2 annotations
    dt_isdb_results = pd.concat([dt_isdb_results, df_MS1])
    dt_isdb_results.rename(columns={'msms_score': 'score_input'}, inplace=True)

    print('Number of annotated features after MS1: ' + str(len(df_MS1['feature_id'].unique())))
    print('Total number of MS1 and MSMS annotations: ' + str(len(dt_isdb_results)))

    # Rank annotations based on the spectral score
    dt_isdb_results["score_input"] = pd.to_numeric(
        dt_isdb_results["score_input"], downcast="float")
    dt_isdb_results['rank_spec'] = dt_isdb_results[['feature_id', 'score_input']].groupby(
        'feature_id')['score_input'].rank(method='dense', ascending=False)

    dt_isdb_results.reset_index(inplace=True, drop=True)

    # now we merge with the Occurences DB metadata after selection of our columns of interest
    cols_to_use = ['structure_smiles_2D', 'structure_molecular_formula',
                'structure_exact_mass', 'short_inchikey', 'structure_taxonomy_npclassifier_01pathway', 
                'structure_taxonomy_npclassifier_02superclass', 'structure_taxonomy_npclassifier_03class',
                'organism_name', 'organism_taxonomy_ottid', 'organism_taxonomy_01domain', 'organism_taxonomy_02kingdom', 'organism_taxonomy_03phylum',
                'organism_taxonomy_04class', 'organism_taxonomy_05order', 'organism_taxonomy_06family', 'organism_taxonomy_08genus',
                'organism_taxonomy_09species', 'organism_taxonomy_10varietas' ]

    dt_isdb_results = pd.merge(
        left=dt_isdb_results, right=db_metadata[cols_to_use], left_on='short_inchikey', right_on='short_inchikey', how='outer')
    dt_isdb_results.dropna(subset=['feature_id'], inplace=True)

    print('''
    Taxonomically informed reponderation
    ''')
    
    # Check if sample has a valid biosource and if not, harmonize ouput 
    if taxo_metadata['matched_name'][0] == 'None':
        dt_isdb_results['score_taxo'] = 0
        dt_isdb_results["score_input"] = pd.to_numeric(dt_isdb_results["score_input"], downcast="float")
        dt_isdb_results['score_input_taxo'] = dt_isdb_results['score_taxo'] +  dt_isdb_results['score_input']
        dt_isdb_results['rank_spec_taxo'] = dt_isdb_results.groupby(
            'feature_id')['score_input_taxo'].rank(method='dense', ascending=False)
        dt_isdb_results = dt_isdb_results.groupby(["feature_id"]).apply(
            lambda x: x.sort_values(["rank_spec_taxo"], ascending=True)).reset_index(drop=True)
        
    # If valid taxonomy is present for sample, proceed to taxonomical reweighting
    else:
        cols_att = ['query_otol_domain', 'query_otol_kingdom', 'query_otol_phylum', 'query_otol_class',
                'query_otol_order', 'query_otol_family', 'query_otol_genus', 'query_otol_species']
        for col in cols_att:
            dt_isdb_results[col] = taxo_metadata[col][0]
        dt_isdb_results = taxonomical_reponderator(dt_isdb_results, min_score_taxo_ms1)
       
    print('''
    Taxonomically informed reponderation done
    ''')
    
    print('''
    Chemically informed reponderation
    ''')

    dt_isdb_results_chem_rew = chemical_reponderator(clusterinfo_summary, dt_isdb_results, top_N_chemical_consistency)
    
    print('''
    Chemically informed reponderation done
    ''')
    
    # Select only the top N annotations and formatting
    dt_taxo_chemo_reweighed_topN = top_N_slicer(input_df=dt_isdb_results_chem_rew,
                                            top_to_output=top_to_output)   
    df_flat, df_for_cyto = annotation_table_formatter(dt_taxo_chemo_reweighed_topN, min_score_taxo_ms1, min_score_chemo_ms1)

    # Export
    df_flat.to_csv(repond_table_flat_path, sep='\t')
    df_for_cyto.to_csv(repond_table_path, sep='\t')

    #Plotting
    feature_intensity_table_formatted = feature_intensity_table_formatter(feature_table)
    plotter_count(df_flat, treemap_chemo_counted_results_path)
    plotter_intensity(df_flat, feature_intensity_table_formatted, treemap_chemo_intensity_results_path)
    
    # Save params 
    shutil.copyfile(r'configs/user/user.yaml', isdb_config_path)
    shutil.copyfile(r'configs/user/user.yaml', mn_config_path)
    
    print('''
    Finished file: ''' + sample_dir
    )