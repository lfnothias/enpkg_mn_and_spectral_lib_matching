# enpkg_mn_and_spectral_lib_matching
Repository gathering scripts to perform molecular networking and spectral library matching on individual files originating from the Feature-Based Molecular Networking worfklow.  

⚙️ Workflow part of [enpkg_workflow](https://github.com/enpkg/enpkg_workflow).

## Required starting architecture and inputs

```
data/
└─── sample_a_folder/
|	└─── pos/
|	     └─── sample_a_features_quant_pos.csv
|	     └─── sample_a_features_ms2_pos.mgf  
|
└─── sample_b_folder/
|
└─── sample_n_folder/
```
OR (auto-restructure)

```
data/
└─── sample_a_folder/
|     └─── sample_a_features_quant_pos.csv
|     └─── sample_a_features_ms2_pos.mgf   
|
└─── sample_b_folder/
|
└─── sample_n_folder/
```


## 1. Clone repository and install environment

1. Clone this repository.
2. Create environment: 
```
conda env create -f environment.yml
```
3. Activate environment:  
```
conda activate enpkg_mn_spectral_lib_matching
```

## 2. Get a spectral library
Download a spectral library of fragmentation spectra (positive and/or negative mode) in .msp format.
For example with GNPS:

```console
wget https://external.gnps2.org/gnpslibrary/ALL_GNPS.msp
```

Move the spectral library into:
```
mv ALL_GNPS.msp db_spectra/
```

## 3. Adapt parameters and launch the process! 🚀

1. Copy and rename the parameters file <code>configs/default/default.yaml</code> into <code>configs/user/user.yaml</code>
2. Modifiy the ``user.yaml`` file according to your needs (especially the paths).
3. Launch the process:

#### Running it
```
python src/enpkg_mn_and_matching_parrallel.py -p data_pos -ion pos -c 10 -r
```

```
python src/enpkg_mn_and_matching_parrallel.py -p data_neg -ion neg -c 10 -r
```

#### Parameters
```
'-p', '--sample_dir_path', required=True, help='The path to the directory where samples folders to process are located'
```
```
'-ion', '--ionization_mode', required=True, choices=['pos', 'neg'], help='The ionization mode to perform spectral library matching'
```
```
'-c', '--cpus', Number of cpu to use. Default is 80% of available CPUs.'
```
```
'-r', '--recompute', Recompute even if the files are already present
```                     
##  Target architecture

```
data/
└─── sample_a_folder/
|     └───  pos/
|			  └───  molecular_network/
|	            └───  sample_a_folder_mn.graphml
|	            └───  sample_a_folder_mn.metadata
				  └───  config.yaml
|		     └───  spectral_lib_matching/
|		         └───  sample_a_folder_lib_results_pos.tsv
|		         └───  sample_a_folder_lib_results_final_pos.tsv
|		     └───  sample_a_features_quant_pos.csv
|		     └───  sample_a_features_ms2_pos.csv 
|
└─── sample_b_folder/
|
└─── sample_n_folder/
```
