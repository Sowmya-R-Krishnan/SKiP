#Python script to write the mol files using RDKit ETKDG method of generating 3D coordinates

import os
import sys
import csv
import pandas as pd

from rdkit import Chem
from rdkit.Chem import AllChem

smilesfile = sys.argv[1]
savepath = sys.argv[2]

#XXX: python 1_generate_mol_files.py "./data/kcat_independent_test_set_v1.csv" "./Sample_output/Ligands/"

df = pd.read_csv(smilesfile, sep="\t", header=0)

for i, row in df.iterrows():
	try:
		smiles = row["Substrate_SMILES"]
		mol_id = row["Mol_ID"]
		rdkitmol = Chem.MolFromSmiles(smiles)
		mol2 = Chem.AddHs(rdkitmol)
		AllChem.EmbedMolecule(mol2, randomSeed=0xf00d)  #Random seed is for reproducibility
		AllChem.MMFFOptimizeMolecule(mol2)

		fname = savepath+str(mol_id)+".sdf"
		writer = Chem.SDWriter(fname)
		writer.write(mol2)
	except:
		fname = savepath+str(mol_id)+".sdf"
		print(smiles+"\t"+str(mol_id))
		os.system("obabel -:\""+smiles+"\" -osdf -O "+fname+" --gen3d -h")
		continue

		

































