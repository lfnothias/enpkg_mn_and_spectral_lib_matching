from matchms.importing import load_from_msp
from matchms.importing import load_from_mgf
from matchms.filtering import default_filters
from matchms.exporting import save_as_mgf
from matchms.exporting import save_as_msp

def load_spectral_db(path_to_db):
    """Load and clean metadata from a .msp spectral database

    Args:
        path_to_db (str): Path to the .msp file

    Returns:
        list: List of matchms spectra object
    """    
    
    print('''
    Cleaning the spectral database metadata fields
    ''')  
    spectrums_db = list(load_from_msp(path_to_db))
    spectrums_db = [default_filters(s) for s in spectrums_db]
    print(f'''
    A total of {len(spectrums_db)} clean spectra were found in the spectral library
    ''')
    return spectrums_db

def load_clean_spectral_db(path_to_db):
    """Load and clean metadata from a .msp spectral database

    Args:
        path_to_db (str): Path to the .msp file

    Returns:
        list: List of matchms spectra object
    """    
    
    print('''
    Loading the spectral library ... this can take minutes but will be done once at the beginning of the process ...
    ''')  
    spectrums_db = list(load_from_msp(path_to_db))
    print(f'''
    A total of {len(spectrums_db)} clean spectra were found in the spectral library
    ''')
    return spectrums_db


def save_spectral_db(spectrums_db, output_path):
    """Save a clean spectral db as .msp from a matchms object

    Args:
        spectrums_db (var): cleaned spectrums_db (e.g. output of load_spectral_db())
        output_path (str): path to output

    Returns:
        list: List of matchms spectra object
    """    
    
    print('''
    Saving the spectral database (it can take minutes ...)
    ''')  
    save_as_msp(spectrums_db, output_path)

    print(f'''
    A total of {len(spectrums_db)} clean spectra were found in the spectral library and saved as {output_path}
    ''')