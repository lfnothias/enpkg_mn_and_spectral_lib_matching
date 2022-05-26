import os
import sys
import glob

def remove_prefix(text, prefix):
    return text[(len(prefix)+1):] if text.startswith(prefix) else text

def main():

    samples_dir = [x[0] for x in os.walk(os.getcwd())]
    root_dir = os.getcwd()
    samples_dir.remove(os.getcwd())

    for sample_dir in samples_dir:

        os.chdir(sample_dir)
        sample_name = remove_prefix(sample_dir, root_dir)
        mgf_file = glob.glob('*.mgf')[0]
        new_mgf_file = sample_name + '_features_ms2_pos.mgf'
        os.rename(mgf_file, new_mgf_file)

        quant_file = glob.glob('*.csv')[0]
        new_quant_file = sample_name + '_features_quant_pos.csv'
        os.rename(quant_file, new_quant_file)

        memo_file = glob.glob('*.csv')[1]
        new_memo_file = sample_name + '_memo_pos.csv'
        os.rename(memo_file, new_memo_file)
              
if __name__ == "__main__":
    main()