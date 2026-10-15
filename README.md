# SKiP: A Structure-Oriented Kinetics Prediction Suite for Enzyme-Substrate Complexes
<b>Authors: Sowmya Ramaswamy Krishnan*, Nishtha Pandey, Arijit Roy*</b>

Automation scripts used for construction of the SKiP suite are arranged into a pipeline. A folder for model training and another folder for inference from pre-trained models is provided. Sample commands to run each script in the preferred order of execution is also provided in a README file within each folder. Models in the SKiP suite were trained using our SKiD dataset (https://doi.org/10.5281/zenodo.15355030). Please contact the corresponding authors in case of any queries or concerns.

<h3>Requirements - Preferably a conda/pip virtual environment with all these packages installed</h3>
<ul>
<li>python>=3.10.13</li>
<li>pandas>=2.3.3</li>
<li>numpy>=1.25.0</li>
<li>scipy>=1.12.0</li>
<li>scikit-learn>=1.7.2</li>
<li>xgboost>=3.0.0</li>
<li>lightgbm>=4.6.0</li>
<li>pickle</li>
<li>matplotlib</li>
<li>seaborn</li>
<li>rdkit-pypi>=2022.9.5</li>
<li>openbabel==3.1.1</li>
<li>biopython>=1.83</li>
<li>openmm>=8.1.1</li>
<li>pdbfixer>=1.9</li>
<li>mordred>=1.2.0</li>
<li>mmseqs2>=16.747c6</li>
<li>mapie>=1.2.0</li>
</ul>

# Other additional stand-alone programs required
* GNINA (pre-compiled binary) - https://github.com/gnina/gnina/releases/tag/v1.3
* FASPR (CLI compiled from source) - https://github.com/tommyhuangthu/FASPR
* PDB2PQR (CLI compiled from source) - https://github.com/Electrostatics/pdb2pqr
* PyMOL (Python bindings) - https://pymol.org/
* NACCESS (CLI compiled from source) - https://www.bioinf.manchester.ac.uk/naccess/

# Usage disclaimer
This is a minimal version of the SKiP model development and inference source code necessary to reproduce the results shown in the manuscript. Wherever possible, appropriate sample input and output files are provided for user reference. Any changes made to the source code (except paths to stand-alone programs) are done at your own risk. The authors will not be liable to any discrepancies observed in the results due to changes made to the source code.

# Data
Sample input and output files and any data used from external databases (PDB, UniProtKB, EMBL Cofactor database etc.,) are provided under the `data` folder of each sub-folder. The PDB files of all enzyme-substrate complexes used during for feature calculations are not included and must be downloaded from the SKiD Zenodo repository prior to running the scripts.

# Code usage
Detailed instructions on how to use the codes are provided under each sub-folder of the repository. For any queries related to code usage, contact the corresponding author for more information.

# Copyright Notice
SKiP code repository is a TCS proprietary resource and should be used for academic purposes only. The contents of this repository should not be used for any commercial purpose without the consent of ALL the authors involved. By downloading and utilizing the scripts, the user consents that any and all Intellectual Property derived from the SKiP code repository is fully owned by TCS in the associated jurisdictions. SKiP code repository usage without citation will be considered illegal.

# Contact Us
For further queries related to code usage, please write to us: roy.arijit3@tcs.com, sowmya.rk1@tcs.com

# Citation
Please cite this article if you use the codes in this repository for your research: 

# License: Creative Commons Attribution Non Commercial No Derivatives 4.0 International











