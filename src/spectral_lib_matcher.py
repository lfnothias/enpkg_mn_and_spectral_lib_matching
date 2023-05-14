import os
import pandas as pd
from tqdm.contrib import tzip
from matchms.filtering import default_filters
from matchms.filtering import normalize_intensities
from matchms.filtering import select_by_intensity
from matchms.filtering import select_by_mz
from matchms.similarity import PrecursorMzMatch
from matchms import calculate_scores
from matchms.similarity import CosineGreedy
from matchms.logging_functions import set_matchms_logger_level
import contextlib
import os
import sys

# See https://github.com/matchms/matchms/pull/271
set_matchms_logger_level("ERROR")

def peak_processing(spectrum):
    spectrum = default_filters(spectrum)
    spectrum = normalize_intensities(spectrum)
    spectrum = select_by_intensity(spectrum, intensity_from=0.01)
    spectrum = select_by_mz(spectrum, mz_from=10, mz_to=1000)
    return spectrum


def spectral_matching(spectrums_query, db_clean, parent_mz_tol,
                      msms_mz_tol, min_cos, min_peaks, output_file_path):
    """Performs spectra matching between query spectra and a database using cosine score

    Args:
        spectrums_query (list): List of matchms spectra objects to query
        db_clean (list): List of reference matchms spectra objects 
        parent_mz_tol (float): Precursor m/z tolerance in Da for matching
        msms_mz_tol (float): m/z tolerance in Da for matching fragments
        min_cos (float): minimal cosine score
        min_peaks (int): minimum number of matching fragments
        output_file_path (str): path to write results
    """
    if os.path.exists(output_file_path):
        os.remove(output_file_path)

    spectrums_query = [peak_processing(s) for s in spectrums_query]

    similarity_score = PrecursorMzMatch(tolerance=parent_mz_tol, tolerance_type="Dalton")
    chunks_query = [spectrums_query[x:x+1000] for x in range(0, len(spectrums_query), 1000)]

    data = []
    for chunk in chunks_query:
        scores = calculate_scores(chunk, db_clean, similarity_score)
        idx_row = scores.scores[:, :][0]
        idx_col = scores.scores[:, :][1]
        import numpy as np

        scans_id_map = {i: int(s.metadata['scans']) for i, s in enumerate(chunk)}

        cosinegreedy = CosineGreedy(tolerance=msms_mz_tol)
        
        for x, y in zip(idx_row, idx_col):
            if x < y:
                msms_score, n_matches = cosinegreedy.pair(chunk[x], db_clean[y])[()]

                if msms_score > min_cos and n_matches > min_peaks:
                    feature_id = scans_id_map[x]
                    
                    data.append({'msms_score': msms_score,
                                 'matched_peaks': n_matches,
                                 'feature_id': feature_id,
                                 'reference_id': y + 1,
                                 'adduct': db_clean[y].get("precusortype"),
                                 'charge': db_clean[y].get("charge"),
                                 'ionmode': db_clean[y].get("ionmode"),
                                 'instrument': db_clean[y].get("instrument"),
                                 'instrument_type': db_clean[y].get("instrumenttype"),  
                                 'comment': db_clean[y].get("comment"),  
                                 'inchikey': db_clean[y].get("inchikey"),
                                 'inchi': (db_clean[y].get("inchi") or '').lstrip('\'"').rstrip('\'"'),
                                 'smiles': db_clean[y].get("smiles"),
                                 'compound_name': db_clean[y].get("compound_name")})

    if len(data) > 0:
        df = pd.DataFrame(data)
        df['Spectral_library_ID'] = df['comment'].str.extract(r'DB#=(.*?);')
        df['Spectral_library'] = df['comment'].str.extract(r'origin=(.*?)(?:;|$)')
    else:
        print("             No spectral library annotations.")
        # Create an empty DataFrame with the same structure
        df = pd.DataFrame(columns=['msms_score', 'matched_peaks', 'feature_id', 'reference_id',
                                   'adduct', 'charge', 'ionmode', 'instrument', 'instrument_type', 
                                   'comment', 'inchikey', 'inchi', 'smiles', 'compound_name',
                                   'Spectral_library_ID', 'Spectral_library'])
    os.makedirs(os.path.dirname(output_file_path), exist_ok=True)
    df.to_csv(output_file_path, mode='a', header=not os.path.exists(output_file_path), sep='\t')






