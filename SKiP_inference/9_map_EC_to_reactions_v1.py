#Program to perform substructure matching between reaction SMILES and SKiD substrates

import os
import sys
import csv
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs

infile = sys.argv[1]
rxns = sys.argv[2]
outfile = sys.argv[3]

#XXX: python 9_map_EC_to_reactions_v1.py "./data/kcat_independent_test_set_v1.csv" "../SKiP_training/data/MetaCyc_reaction_SMILES_v2.csv" "./Sample_output/kcat_MetaCyc_rxns_v1.csv"

skid = pd.read_csv(infile, sep="\t", header=0)
rxn_df = pd.read_csv(rxns, sep="\t", header=0)
#print(skid.columns)  #['Entry_ID', 'EC_number', 'Substrate', 'UniProt_ID', 'Protein_file', 'Organism_name', 'Km_value (mM)', 'Mutant', 'Mutation', 'pH', 'Temperature', 'References', 'Site_type', 'Substrate_SMILES', 'Mol_file', 'pKm_value', 'UniProt_URL', 'BRENDA_URL', 'Parent_PDB_structure', 'ES_complex_modelling']

out = open(outfile, "w")
print("Entry_ID\tEC_number\tPDB_ID\tpH\tSubstrate_SMILES\tMol_ID\tMaximum_MetaCyc_match\tReaction_ID\tEC_number\tType\tTC_value", file=out)

done = []
for i, row in skid.iterrows():
	if(row["Entry_ID"] not in done):
		substrate = row["Substrate_SMILES"]
		ec_num = row["EC_number"].split("; ")[0]
		
		ec_parts = ec_num.split(".")
		ec_new = ec_parts[0]+"."+ec_parts[1]+"."+ec_parts[2]
		sub_df = rxn_df[rxn_df["EC_number"].str.startswith("EC-"+ec_new)]
		
		sims = []
		max_tc = ""
		maximum = 0
		if(len(sub_df.index)>0):
			for j, r in sub_df.iterrows():
				fp1 = AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(substrate), 2, nBits=1024)
				reactants = r["Reactants"].split(", ")
				products = r["Products"].split(", ")
				
				for react in reactants:
					fp2 = AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(react), 2, nBits=1024)
					tc = DataStructs.TanimotoSimilarity(fp1, fp2)
					sims.append((react, r["Reaction_ID"], r["EC_number"], tc, "Reactant"))
					if(tc > maximum):
						max_tc = (react, r["Reaction_ID"], r["EC_number"], tc, "Reactant")
						maximum = tc
						
				for prod in products:
					fp2 = AllChem.GetMorganFingerprintAsBitVect(Chem.MolFromSmiles(prod), 2, nBits=1024)
					tc = DataStructs.TanimotoSimilarity(fp1, fp2)
					sims.append((prod, r["Reaction_ID"], tc, "Product"))
					if(tc > maximum):
						max_tc = (prod, r["Reaction_ID"], r["EC_number"], tc, "Product")
						maximum = tc
			
			if(max_tc!=""):
				print(row["Entry_ID"]+"\t"+str(row["EC_number"])+"\t"+str(row["PDB_ID"])+"\t"+str(row["pH"])+"\t"+substrate+"\t"+row["Mol_ID"]+"\t"+max_tc[0]+"\t"+max_tc[1]+"\t"+max_tc[2]+"\t"+max_tc[-1]+"\t"+str(max_tc[3]), file=out)
				done.append(row["Entry_ID"])
			else:
				continue
		else:
			print(row["EC_number"])


























