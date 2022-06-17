from matchms.importing import load_from_mgf
from matchms.filtering import default_filters

def load_spectral_db(path_to_db):
    """Load and clean metadata from a .mgf spectral database

    Args:
        path_to_db (str): Path to the .mgf file

    Returns:
        list: List of matchms spectra object
    """    
    
    print('''
    Cleaning the spectral database metadata fields
    ''')  
    spectrums_db = list(load_from_mgf(path_to_db))
    spectrums_db = [default_filters(s) for s in spectrums_db]
    print(f'''
    A total of {len(spectrums_db)} clean spectra were found in the spectral library
    ''')
    return spectrums_db
