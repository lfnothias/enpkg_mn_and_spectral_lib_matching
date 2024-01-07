from pathlib import Path
import pandas as pd
import shutil
import sys
import os
from tqdm import tqdm

# We check for pos or neg presence and we construct from the polarity argument

def check_for_ionmode_folder_and_restruct_if_needed(path, polarity):
    path = os.path.normpath(path)
    samples_dir = [directory for directory in os.listdir(path) if not directory.startswith('.DS_Store')]

    for directory in tqdm(samples_dir):
        pos_folder = os.path.join(path, directory, 'pos')
        neg_folder = os.path.join(path, directory, 'neg')

        if os.path.exists(pos_folder):
            ionization_mode = 'pos'
        elif os.path.exists(neg_folder):
            ionization_mode = 'neg'
        else:
            target_dir = os.path.join(path, directory, polarity)
            if not os.path.exists(target_dir):
                os.makedirs(target_dir)

            dir_path = os.path.join(path, directory)
            for entry in os.scandir(dir_path):
                if entry.is_file(): #and entry.name != 'metadata.tsv':
                    shutil.copy(entry.path, os.path.join(target_dir, entry.name))
            print(target_dir)

# We check for pos or neg presence and we construct from the SIRIUS folder
def check_for_ionmode_folder_and_autorestruct_from_sirius_if_needed(path):
    path = os.path.normpath(path)
    samples_dir = [directory for directory in os.listdir(path) if not directory.startswith('.DS_Store')]
    pos_folder = os.path.join(path, directory, 'pos')
    neg_folder = os.path.join(path, directory, 'neg')

    for directory in tqdm(samples_dir):
        if os.path.exists(pos_folder):
            ionization_mode = 'pos'
        elif os.path.exists(neg_folder):
            ionization_mode = 'neg'
        else:
            csi_path = os.path.join(path, directory, 'compound_identifications.tsv')
            csi_annotations = pd.read_csv(csi_path, sep='\t')

            adduct_counts = csi_annotations['adduct'].value_counts()
            most_frequent_adduct = adduct_counts.idxmax()

            if ']+' in most_frequent_adduct:
                ionization_mode = 'pos'
                print('Auto gave: '+ ionization_mode)
            elif ']-' in most_frequent_adduct:
                ionization_mode = 'neg'
            else:
                raise ValueError('Cannot deduce polarity from the most frequent adduct.')

        target_dir = os.path.join(path, directory, ionization_mode)
        if not os.path.exists(target_dir):
            os.makedirs(target_dir)

        dir_path = os.path.join(path, directory)
        for entry in os.scandir(dir_path):
            if entry.is_file() and entry.name != 'metadata.tsv':
                shutil.move(entry.path, os.path.join(target_dir, entry.name))
        print(target_dir)