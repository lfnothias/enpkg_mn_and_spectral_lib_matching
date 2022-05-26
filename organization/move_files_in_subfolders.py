import os
import shutil

def main():

    cwd = os.getcwd()

    for root, subFolders, files in os.walk(cwd):
        for file in files:
            if file.endswith('_sirius_neg.mgf'):
                subFolder = os.path.join(cwd, file[:10])
                shutil.move(os.path.join(root, file), subFolder)
            elif file.endswith('_ms2_neg.mgf'):
                subFolder = os.path.join(cwd, file[:10])
                shutil.move(os.path.join(root, file), subFolder)
            elif file.endswith('_quant_neg.csv'):
                subFolder = os.path.join(cwd, file[:10])
                shutil.move(os.path.join(root, file), subFolder)
              
if __name__ == "__main__":
    main()