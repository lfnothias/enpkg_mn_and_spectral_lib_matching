# required libraries

import pandas as pd
import numpy as np
import zipfile
import glob
import os
import sys
import shlex
import subprocess
import json
from pandas import json_normalize
from tqdm import tqdm
#from pivottablejs import pivot_ui

def cluster_counter(clusterinfo_summary_file):
    """ Count the numbers of nodes per component index in a molecular network

    Args:
        clusterinfo_summary_file (dataframe) : a molecular network clusterinfo_summary_file
    Returns:
        cluster_count (dataframe): a dataframe with the number of nodes per component index
    """


    cluster_count = clusterinfo_summary_file.drop_duplicates(
        subset=['feature_id', 'component_id']).groupby("component_id").count()
    cluster_count = cluster_count[['feature_id']].rename(
        columns={'feature_id': 'ci_count'}).reset_index()
    return cluster_count


def top_N_slicer(dt_isdb_results, top_to_output):

    """ Keeps only the top N candidates out of an annotation table and sorts them by rank

    Args:
        dt_isdb_results (dataframe) : a annotation table
        top_to_output (integer): Top N of candidate to keep
    Returns:
        dt_isdb_results_chem_rew (dataframe): a dataframe with the top N annotation ordered by final rank
    """

    dt_isdb_results_chem_rew = dt_isdb_results.loc[(
        dt_isdb_results.rank_final <= int(top_to_output))]
    dt_isdb_results_chem_rew[["feature_id", "rank_final", "component_id"]] = dt_isdb_results_chem_rew[[
        "feature_id", "rank_final", "component_id"]].apply(pd.to_numeric, downcast='signed', axis=1)
    dt_isdb_results_chem_rew = dt_isdb_results_chem_rew.sort_values(
        ["feature_id", "rank_final"], ascending=(False, True))

    return dt_isdb_results_chem_rew


def annotation_table_formatter(dt_input, min_score_taxo_ms1, min_score_chemo_ms1):
    """ A bunche of formatter frunctions for output

    Args:
        dt_input (dataframe) : an annotation table
        keep_lowest_taxon (bool): wether to export only the lowest matched taxon or all the taxonomy for annotations
        top_to_output (integer): Top N of candidate to keep
    Returns:
        dt_isdb_results_chem_rew (dataframe): a dataframe with the top N annotation ordered by final rank
    """

    # Here we would like to filter results when short IK are repeated for the same feature_id at the same final rank

    dt_input = dt_input.drop_duplicates(
        subset=['feature_id', 'short_inchikey'], keep='first')

    dt_input = dt_input.astype(
        {'feature_id': 'int64'})

    dt_input['lowest_matched_taxon'] = dt_input['matched_species']
    dt_input['lowest_matched_taxon'] = dt_input['lowest_matched_taxon'].replace(
        'nan', np.NaN)
    col_matched = ['matched_genus', 'matched_family', 'matched_order',
                    'matched_order', 'matched_phylum', 'matched_kingdom', 'matched_domain']
    for col in col_matched:
        dt_input[col] = dt_input[col].replace(
            'nan', np.NaN)
        dt_input['lowest_matched_taxon'].fillna(
            dt_input[col], inplace=True)

    annot_attr = ['rank_spec', 'score_input', 'libname', 'short_inchikey', 'structure_smiles_2D', 'structure_molecular_formula', 'adduct',
                    'structure_exact_mass', 'structure_taxonomy_npclassifier_01pathway', 'structure_taxonomy_npclassifier_02superclass',
                    'structure_taxonomy_npclassifier_03class',
                    'query_otol_species', 'lowest_matched_taxon', 'score_taxo', 'score_max_consistency', 'final_score', 'rank_final']

    comp_attr = ['component_id', 'structure_taxonomy_npclassifier_01pathway_consensus', 'freq_structure_taxonomy_npclassifier_01pathway',
                 'structure_taxonomy_npclassifier_02superclass_consensus',
                 'freq_structure_taxonomy_npclassifier_02superclass', 'structure_taxonomy_npclassifier_03class_consensus', 'freq_structure_taxonomy_npclassifier_03class']

    col_to_keep = ['feature_id'] + comp_attr + annot_attr

    # We add the min chemo score at this step
    dt_input = dt_input[
        ((dt_input['score_taxo'] >= min_score_taxo_ms1) & (dt_input['score_max_consistency'] >= min_score_chemo_ms1)) | (
            dt_input['libname'] == 'ISDB')]
    dt_output_flat = dt_input[col_to_keep]

    # Cytoscape formatting 

    all_columns = list(dt_input) # Creates list of all column headers   
    dt_input[all_columns] = dt_input[all_columns].astype(str)
    gb_spec = {c: '|'.join for c in annot_attr}

    for c in comp_attr:
        gb_spec[c] = 'first'

    dt_output_cyto = dt_input.groupby('feature_id').agg(gb_spec)
    dt_output_cyto.reset_index(inplace=True)

    return dt_output_flat, dt_output_cyto


