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
```console
mv ALL_GNPS.msp db_spectra/
```

## 3. Adapt parameters and launch the process! 🚀

1. Copy and rename the parameters file <code>configs/default/default.yaml</code> into <code>configs/user/user.yaml</code>
2. Modifiy the ``user.yaml`` file according to your needs (especially the paths).
3. Launch the process:
```
python src/enpkg_mn_and_matching.py -p data -ion pos
```

##  Target architecture

```
data/
└─── sample_a_folder/
|     └───  pos/
|			  └───  molecular_network/
|	            └───  sample_a_folder_mn.graphml
|	            └───  sample_a_folder_mn.metadata
|		     └───  spectral_lib_matching/
|		         └───  sample_a_folder_db_results.tsv
|		         └───  sample_a_folder_db_results_final.tsv
|		     └───  sample_a_features_quant_pos.csv
|		     └───  sample_a_features_ms2_pos.csv 
|
└─── sample_b_folder/
|
└─── sample_n_folder/
```
