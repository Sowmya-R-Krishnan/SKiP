# SKiP: A Structure-Oriented Kinetics Prediction Suite for Enzyme-Substrate Complexes
<b>Authors: Sowmya Ramaswamy Krishnan*, Nishtha Pandey, Arijit Roy*</b>

Sample commands for running all the scripts in this sub-folder are provided with kcat dataset as the example. Users should replicate the same commands for Km dataset for building the model. Please contact the corresponding authors in case of any queries or concerns.

# Scripts for training the machine learning models in SKiP suite
1. compute_model_features_v1.py - Script to compute global enzyme features, local enzyme-substrate interaction features and substrate features for a set of enzyme-substrate complex structures.
2. map_EC_to_reactions_v1.py - Script to map EC numbers to the dataset of MetaCyc reactions for RDM resolution.
3. get_reaction_frequency_v1.py - Script to use Tanimoto coefficient cut-off between substrate and MetaCyc reactants to map the top 3 RDMs corresponding to the EC-substrate pair.
4. compute_R_part_fingerprint_v1.py - Script to compute the RDM fingerprints for the input dataset.
5. combine_feat_files_v1.py - Script to combine global features, local features, substrate features, RDM features and the predictive variable (pkcat or pKm) to form the final feature vectors for model development.
6. split_dataset_v1.py - Script to perform random split or EC-wise split of the dataset into training and test sets.
7. train_model_v1.py - Script to train various machine learning models for kinetic parameter prediction in SKiP.
8. predict_uncertainty_v1.py - Script to delineate the aleatoric and epistemic uncertainty contributions from the pre-trained ML model for each data point in a given dataset.
9. combine_uncertainty_preds_v1.py - Script to combine aleatoric and epistemic uncertainty estimates into a single uncertainty value and its standard deviation.

# Requirements - Preferably a conda environment with all these packages installed
* python>=3.10.13
* pandas>=2.3.3
* numpy>=1.25.0
* scipy>=1.12.0
* scikit-learn>=1.7.2
* xgboost>=3.0.0
* lightgbm>=4.6.0
* pickle
* matplotlib
* seaborn
* rdkit-pypi>=2022.9.5
* openbabel==3.1.1
* biopython>=1.83
* openmm>=8.1.1
* pdbfixer>=1.9
* mordred>=1.2.0
* mmseqs2>=16.747c6
* mapie>=1.2.0

# Other additional stand-alone programs required
* PyMOL (Python bindings) - https://pymol.org/
* NACCESS (CLI compiled from source) - https://www.bioinf.manchester.ac.uk/naccess/

# Data
Sample input and output files and any data used from external databases (PDB, UniProtKB, EMBL Cofactor database etc.,) are provided under `data`, `Sample_structures` and `Sample_output` folders. The pickle files of final models used in SKiP study, including uncertainty prediction results from the MAPIE package, are provided under the `Final_SKiP_models` folder. The PDB files of all the enzyme-substrate complexes used during the calculations are not included and must be downloaded from the SKiD Zenodo repository prior to running the scripts.

# Code usage - Sample commands
```
* python 1_compute_model_features_v1.py "./data/kcat_with_Km_WT_BS_data_logscale_final_v1.csv" "./Sample_structures/SKiD_kcat_archive/" "./Sample_structures/SKiD_kcat_archive/" "./Sample_output/Sample_local_feats_v1.csv" "./Sample_output/Sample_global_feats_v1.csv" "./Sample_output/Sample_substrate_feats_v1.csv" "./" "./data/Pathway_mapping_final.csv" <naccess binary path>
* python 2_map_EC_to_reactions_v1.py "./data/kcat_with_Km_WT_BS_data_logscale_final_v1.csv" "./data/MetaCyc_reaction_SMILES_v2.csv" "./Sample_output/kcat_dataset_MetaCyc_rxns_sample.csv"
* python 3_get_reaction_frequency_v1.py "./Sample_output/kcat_dataset_MetaCyc_rxns_sample.csv" "./data/uniqueRDM_db.csv" "./data/cofactor_EC_v02 1.csv" "./Sample_output/kcat_dataset_RDMs_unmapped_sample.csv" "./Sample_output/kcat_dataset_RDMs_mapped_sample.csv"
* python 4_compute_R_part_fingerprint_v1.py "./data/kcat_with_Km_WT_BS_data_logscale_final_v1.csv" "./Sample_output/kcat_dataset_RDMs_mapped_sample.csv" "./Sample_output/kcat_dataset_RDMs_unmapped_sample.csv" "./data/Atom_types.xlsx" "./Sample_output/kcat_sample_RDM_fingerprints_v1.csv" "./data/In_house_RDMs.csv"
* python 5_combine_feat_files_v1.py "./data/kcat_header.csv" "./data/kcat_with_Km_WT_BS_data_logscale_final_v1.csv" "./Sample_output/Sample_kcat_local_feats_v1.csv" "./Sample_output/Sample_kcat_global_feats_v1.csv" "./Sample_output/Sample_kcat_substrate_feats_v1.csv" "./Sample_output/kcat_sample_RDM_fingerprints_v1.csv" "./Sample_output/kcat_sample_final_feats_v1.csv"
* python 6_split_dataset_v1.py "./Sample_output/kcat_sample_final_feats_v1.csv" "random" "kcat" "./data/kcat_with_Km_WT_BS_data_logscale_final_v1.csv" "./Sample_output/Random_split/"
* python 6_split_dataset_v1.py "./Sample_output/kcat_sample_final_feats_v1.csv" "strat_ec_1" "kcat" "./data/kcat_with_Km_WT_BS_data_logscale_final_v1.csv" "./Sample_output/EC_wise_split/"
* python 7_train_model_v1.py "./Sample_output/Random_split/" "kcat" 3 "./Sample_output/Random_split/kcat_LGBM_v1.pkl"
* python 8_predict_uncertainty_v1.py "./Sample_output/Random_split/kcat_X_test.csv" 3 "kcat" "./Final_SKiP_models/kcat_LGBM_10F_CV/" "./Final_SKiP_models/kcat_LGBM_MAPIE_results_v1.pkl" "./Sample_output/kcat_X_test_uncertainty_sample.csv" "./Sample_output/Random_split/kcat_Y_test.csv"
* python 9_combine_uncertainty_preds_v1.py "./Sample_output/kcat_X_test_uncertainty_sample.csv" "./Sample_output/kcat_X_test_uncertainty_combined_sample.csv"
```

# Miscellaneous
All inputs to the codes are provided as command-line arguments. Any output folder required to write the results must be created by the user before running the scripts. Only paths to stand-alone programs like NACCESS or PyMOL python bindings should be modified in the code to ensure reproducibility. Modification of any other parameters is at your own risk.

* Path to the NACCESS binary must be added to the system path to run Script 1.
* Entry_ID is an important column in the input files and the structure archives must be built based on Entry_ID to label folders containing structure information. Please verify the correctness of input files against the datasets provided in the `data` folder.
* TODO statements are used to highlight lines in the code which must be commented/uncommented for kcat and Km datasets. Especially, scripts 5 and 7 need to be carefully edited depending on the input and model requirements.

# License: Creative Commons Attribution Non Commercial No Derivatives 4.0 International










