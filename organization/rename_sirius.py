import os
import sys
import glob

def main():
    
    files = glob.glob("*.mgf")
    for file in files:        
        new_mgf_file = file[0:10] + '_sirius_pos.mgf'
        os.rename(file, new_mgf_file)
              
if __name__ == "__main__":
    main()