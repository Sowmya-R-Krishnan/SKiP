# SKiP: A Structure-Oriented Kinetics Prediction Suite for Enzyme-Substrate Complexes
<b>Authors: Sowmya Ramaswamy Krishnan*, Nishtha Pandey, Arijit Roy*</b>

Before running the scripts inside this folder, kindly prepare the input file in the same format as given in `data` folder (TSV format). The scripts assume that the SIFTS XML files (*.xml), input PDB files (*.pdb) and the reference ligands for docking (*_orig.sdf) calculations are already available with the user in the `data` folder. Ensure that these files are available before running the scripts to avoid errors. If you are unsure how to prepare the structure files, please refer to the SKiD pipeline (https://github.com/Sowmya-R-Krishnan/SKiD). Sample commands for running all the scripts in this sub-folder are provided with kcat dataset as the example. Users should replicate the same commands for Km dataset for building the model. Please contact the corresponding authors in case of any queries or concerns.

# Scripts for running inference using the pre-trained machine learning models in SKiP suite
1. generate_mol_files.py - Script to generate 3D structures of substrates in SDF format using RDKit or OpenBabel CLI, in case RDKit fails.
2. map_brenda_mut_to_bs_v1.py - Script to map mutations listed in BRENDA to their corresponding positions in PDB using EMBL SIFTS mapping.
3. prepare_FASPR_input_v1.py - Script to prepare input file for FASPR mutant modelling calculations.
4. FASPR_automation_v1.py - Script to run the FASPR tool to model enzyme mutants for docking calculations.
5. run_pdbfixer_and_pymol_v1.py - Script to remove any hydrogens present and fill in missing heavy atoms and residues for pH-conditioned protonation.
6. protonate_propka_v1.py - Script to use the Propka program from PDB2PQR package to protonate the enzyme structures based on the input pH.
7. run_docking_calc_v1.py - Script to run docking calculations using GNINA with CNN-rescoring and multi-pose sampling (n=10).
8. compute_model_features_v1.py - Script to compute global enzyme features, local enzyme-substrate interaction features and substrate features for a set of enzyme-substrate complex structures.
9. map_EC_to_reactions_v1.py - Script to map EC numbers to the dataset of MetaCyc reactions for RDM resolution.
10. get_reaction_frequency_v1.py - Script to use Tanimoto coefficient cut-off between substrate and MetaCyc reactants to map the top 3 RDMs corresponding to the EC-substrate pair.
11. compute_R_part_fingerprint_v1.py - Script to compute the RDM fingerprints for the input dataset.
12. create_final_datasets_v1.py - Script to combine global features, local features, substrate features and RDM features to form the final feature vectors for model inference.
13. predict_Km_vals_v1.py - Script to predict pKm values from the pre-trained SKiP model.
14. predict_kcat_vals_v1.py - Script to predict pkcat values from the pre-trained SKiP model.
15. compute_metrics_v1.py - Script to compute model performance metrics if the experimental kinetic parameters are available with the user.

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
* GNINA (pre-compiled binary) - https://github.com/gnina/gnina/releases/tag/v1.3
* FASPR (CLI compiled from source) - https://github.com/tommyhuangthu/FASPR
* PDB2PQR (CLI compiled from source) - https://github.com/Electrostatics/pdb2pqr
* PyMOL (Python bindings) - https://pymol.org/
* NACCESS (CLI compiled from source) - https://www.bioinf.manchester.ac.uk/naccess/

# Data
Sample input and output files and any data used from external databases (PDB, UniProtKB, EMBL Cofactor database, EMBL SIFTS etc.,) are provided under `data` and `Sample_output` folders, respectively. The pickle files of final models used in SKiP study, including uncertainty prediction results from the MAPIE package, are provided under the `Final_SKiP_models` sub-folder under the `SKiP_training` folder. The PDB files, SIFTS XML files and reference ligands used during the calculations must be downloaded and defined by the user prior to running the scripts.

# Code usage - Sample commands
```
* python 1_generate_mol_files.py "./data/kcat_independent_test_set_v1.csv" "./Sample_output/Ligands/"
* python 2_map_brenda_mut_to_bs_v1.py "./data/kcat_independent_test_set_v1.csv" "./data/SIFTS_mapping/" "./Sample_output/kcat_test_set_PDB_muts_v1.csv"
* python 3_prepare_FASPR_input_v1.py "./Sample_output/kcat_test_set_PDB_muts_v1.csv" "./Sample_output/kcat_FASPR_input_v1.csv"
* python 4_FASPR_automation_v1.py "./Sample_output/kcat_FASPR_input_v1.csv" <path to FASPR source code>
* python 5_run_pdbfixer_and_pymol_v1.py "./Sample_output/FASPR_output/" "./Sample_output/Processed_PDBFixer/" "./Sample_output/Processed_PyMOL/"
* python 6_protonate_propka_v1.py "./data/kcat_independent_test_set_v1.csv" "./Sample_output/Processed_PyMOL/" "./Sample_output/Docking_ready/" "./Sample_output/Propka_logs/"
* python 7_run_docking_calc_v1.py "./data/kcat_independent_test_set_v1.csv" "./Sample_output/Docking_ready/" "./Sample_output/Ligands/" "./Sample_output/Docking_results/" <path to GNINA executable file>
* python 8_compute_model_features_v4.py "./data/kcat_independent_test_set_v1.csv" "./Sample_output/Docking_ready/" "./Sample_output/Docking_results/" "./Sample_output/Features/kcat_local_feats_v1.csv" "./Sample_output/Features/kcat_global_feats_v1.csv" "./Sample_output/Features/kcat_substrate_feats_v1.csv" "./" "../SKiP_training/data/Pathway_mapping_final.csv" <path to naccess executable>
* python 9_map_EC_to_reactions_v1.py "./data/kcat_independent_test_set_v1.csv" "../SKiP_training/data/MetaCyc_reaction_SMILES_v2.csv" "./Sample_output/kcat_MetaCyc_rxns_v1.csv"
* python 10_get_reaction_frequency_v1.py "./Sample_output/kcat_MetaCyc_rxns_v1.csv" "../SKiP_training/data/uniqueRDM_db.csv" "../SKiP_training/data/cofactor_EC_v02 1.csv" "./Sample_output/kcat_RDMs_unmapped_v1.csv" "./Sample_output/kcat_RDMs_mapped_v1.csv"
* python 11_compute_R_part_fingerprint_v1.py "./data/kcat_independent_test_set_v1.csv" "./Sample_output/kcat_RDMs_mapped_v1.csv" "./Sample_output/kcat_RDMs_unmapped_v1.csv" "../SKiP_training/data/Atom_types.xlsx" "./Sample_output/Features/kcat_RDM_fingerprints_v1.csv" "../SKiP_training/data/In_house_RDMs.csv"
* python 12_create_final_datasets_v1.py "../SKiP_training/Sample_output/Random_split/kcat_X_test.csv" "./data/kcat_independent_test_set_v1.csv" "./Sample_output/Features/kcat_local_feats_v1.csv" "./Sample_output/Features/kcat_global_feats_v1.csv" "./Sample_output/Features/kcat_substrate_feats_v1.csv" "./Sample_output/Features/kcat_RDM_fingerprints_v1.csv" "./Sample_output/Features/kcat_independent_test_set_feats_final_v1.csv"
* python 13_predict_Km_vals_v1.py "./Sample_output/Features/Km_independent_test_set_final_feats_v1.csv" "../SKiP_training/Final_SKiP_models/Km_RF_v1.pkl" 3 "./Sample_output/SKiP_Km_preds_v1.csv"
* python 14_predict_kcat_vals_v1.py "./Sample_output/Features/kcat_independent_test_set_feats_final_v1.csv" "../SKiP_training/Final_SKiP_models/kcat_LGBM_v1.pkl" 3 "./Sample_output/SKiP_kcat_preds_v1.csv"
* python 15_compute_metrics_v1.py "./data/kcat_independent_test_set_v1.csv" "./Sample_output/SKiP_kcat_preds_v1.csv" "pkcat_value" "Pred_pkcat"
```

# Miscellaneous
All inputs to the codes are provided as command-line arguments. Any output folder required to write the results must be created by the user before running the scripts. Only paths to stand-alone programs like FASPR, Propka, NACCESS or PyMOL python bindings should be modified in the code to ensure reproducibility. Modification of any other parameters is at your own risk.

* Folders required for output storage under the `Sample_output` sub-folder such as, `FASPR_input`, `FASPR_output`, `Ligands` etc must be created by the user before running the appropriate scripts.
* Please verify the correctness of input files against the datasets provided in the `data` folder.
* Path to the source-compiled FASPR tool must be provided as a command line argument to run Script 4.
* Path to the GNINA (v1.1) binary executable must be provided as a command line argument to run Script 7.
* Path to the NACCESS binary must be added to the system path to run Script 8.
* TODO statements are used to highlight lines in the code which must be commented/uncommented for kcat and Km datasets. Especially, script 12 needs to be carefully edited depending on the input requirements.

# License: Creative Commons Attribution Non Commercial No Derivatives 4.0 International










