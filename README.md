# indifiles_annotation
Repository gathering scripts to perform MS/MS annotation using ISDB, MS1 annotation, Molecular Networking and annotations reweighting on a set of individual files.  

▶️ This processing is part of [enpkg_workflow](https://github.com/mandelbrot-project/enpkg_workflow). Data have to be organized using [data_organization](https://github.com/mandelbrot-project/data_organization) and taxonomy resolved using [taxo_enhancer](https://github.com/mandelbrot-project/taxo_enhancer).  

## Required starting architecture

⚠️ Positive mode only for the moment

```
data/
└─── sample_a/
|     └─── sample_a_metadata.tsv
|     └─── pos/
|     |     sample_a_features_quant_pos.csv
|     |     sample_a_features_ms2_pos.mgf
|     |     
|     └─── taxo_output/
|            └─── sample_a_species.json
|            └─── sample_a_taxon_info.json
|            └─── sample_a_taxo_metadata.tsv
|
└─── sample_b/
|
└─── sample_n/
```

## 1. Clone repository and install environment

1. Clone this repository.
2. Create environment: 
```console 
conda env create -f environment.yml
```
3. Activate environment:  
```console 
conda activate indifiles_annotation_env
```

## 2. Get structure-organism pairs and spectral database
Thanks to the [LOTUS initiative](https://lotus.nprod.net/), a large number of NPs with their associated biosources are made available. *In silico* fragmented spectra of these NPs are also available.  
1. Download structure-organism pairs: https://zenodo.org/record/6582124#.YqwzU3ZBxPY
2. Download *in silico* fragmentation spectra: https://zenodo.org/record/5607264#.Yqwwk3ZBxPY
3. Move the structure-organism pairs file into:  
<code>../indifiles_annotation/db_metadata/</code>
3. Move the spectra file into:  
<code>../indifiles_annotation/db_spectra/</code>

## 3. Prepare potential adducts

```console
python src/adducts_formatter.py -p db_metadata/220525_frozen_metadata.csv.gz # Replace according to your version
```
This will create adducts 2 adducts files (pos/neg) and the adducts used (params) in:  
<code>../indifiles_annotation/data_loc/220525_frozen_metadata/</code>

NB: To edit calculated adducts, modify this script according to your needs:  
https://github.com/mandelbrot-project/indifiles_annotation/blob/main/src/adducts_formatter.py

## 3. Adapt parameters and launch the process! 🚀

1. Copy and rename the parameters file <code>../indifiles_annotation/configs/default/default.yaml</code> into <code>../indifiles_annotation/configs/user/user.yaml</code>
2. Modifiy the user.yaml file according to your needs (especially the paths).
3. Launch the process:
```console
python src/nb_indifile.py
```

##  Target architecture

```
data/
└─── sample_a/
|     └───  sample_a_metadata.tsv
|     └─── pos/
|     |     |  sample_a_features_quant_pos.csv
|     |     |  sample_a_features_ms2_pos.mgf
|     |     └─── isdb/
|     |     |      └─── sample_a_isdb_pos.tsv                       # MS2 annootations
|     |     |      └─── sample_a_isdb_reweighted_pos.tsv            # MS2 annotations reweighted, cytoscape ready (1 feature by line)
|     |     |      └─── sample_a_isdb_reweighted_flat_pos.tsv       # MS2 annotations flat (1 annotation by line)
|     |     |      └─── sample_a_treemap_chemo_counted_pos.html     # NPClassifier treemap using annotation count
|     |     |      └─── sample_a_treemap_chemo_intensity_pos.html   # NPClassifier treemap using annotation average intensity
|     |     |      └─── config.yaml                                 # configuration used
|     |     └─── isdb/
|     |     |      └─── sample_a_mn_pos.graphml                     # Moleculat network
|     |     |      └─── sample_a_mn_metadata_pos.tsv                # Moleculat network metadata (component id, precursor m/z)
|     |     |      └─── config.yaml                                 # configuration used
|     └─── taxo_output/
|            └─── sample_a_species.json
|            └─── sample_a_taxon_info.json
|            └─── sample_a_taxo_metadata.tsv
|
└─── sample_b/
|
└─── sample_n/
```
