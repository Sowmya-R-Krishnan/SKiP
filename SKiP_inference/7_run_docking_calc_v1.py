#Program to automate docking of Platinum database ligands with their respective receptors

import os
import sys
import csv
import pandas as pd
import numpy as np

infile = sys.argv[1]
protpath = sys.argv[2]
ligpath = sys.argv[3]
outpath = sys.argv[4]
gnina_path = sys.argv[5]

#XXX: python 7_run_docking_calc_v1.py "./data/kcat_independent_test_set_v1.csv" "./Sample_output/Docking_ready/" "./Sample_output/Ligands/" "./Sample_output/Docking_results/" <path to GNINA executable file>

dock_df = pd.read_csv(infile, sep="\t", header=0)
#print(dock_df.columns)  #['EC_number', 'Substrate', 'UniProt_ID', 'Protein_file', 'Organism_name', 'Km_value', 'Mutant', 'Mutation', 'pH', 'Temperature', 'References', 'Site_type', 'Substrate_SMILES', 'Mol_file']

three_to_one = {"ALA":"A", "GLY":"G", "VAL":"V", "LEU":"L", "ILE":"I", "PHE":"F", "TRP":"W", "TYR":"Y", "HIS":"H", "ARG":"R", "LYS":"K", "ASN":"N", "GLN":"Q", "ASP":"D", "GLU":"E", "SER":"S", "CYS":"C", "THR":"T", "MET":"M", "PRO":"P", "OCS":"C", "CSO":"C", "MSE":"M", "MHS":"H", "CME":"C", "CSD":"C", "TSY":"C", "CSS":"C", "SCY":"C", "PTR":"Y", "DDZ":"A", "JJJ":"C", "MHO":"M", "TPO":"T", "SEP":"S", "KPI":"K", "SLZ":"K", "KPF":"K", "KPY":"K", "74P":"K", "VPV":"K", "KGC":"K", "KYQ":"K", "LYF":"K", "CIR":"R", "SNC":"C", "CSP":"C", "LLP":"K", "KCX":"K"}

for i, row in dock_df.iterrows():
	prot_file = row["PDB_ID"]
	if(row["Mutant"]=="yes"):
		brenda_mut = row["Mutations"].split("/")
		PDB_mut = row["PDB_muts"].split("/")
		
		mut_str = []
		for bmut, pmut in zip(brenda_mut, PDB_mut):  #TODO: Assumes that both mutation lists follow the same order! Make sure its valid!!
			pdb_mut_new = three_to_one[pmut.split("_")[0]]+pmut.split("_")[1]+bmut[-1]
			mut_str.append(row["Chain"]+":"+pdb_mut_new)
			
		prot_file = row["PDB_ID"]+"_"+row["Chain"]+"_"+"+".join(mut_str)
		
	try:
		ph_val = float(row["pH"])
	except:
		ph_val = 7.0

	prot_name = prot_file[0:4]
	prot_file = prot_file+"_"+str(ph_val)+".pdb"
	lig_file = row["Mol_ID"]+".sdf"
	ref_lig = prot_name+"_"+row["Ref_lig"]+"_orig.sdf"
	
	print(prot_file, lig_file, ref_lig)
	
	if(os.path.exists(protpath+prot_file) and os.path.exists(ligpath+lig_file)):
		os.system(gnina_path+" -r "+protpath+prot_file+" -l \""+ligpath+lig_file+"\" --autobox_ligand="+protpath+ref_lig+" -o \""+outpath+prot_file.replace(".pdb","")+"_"+lig_file.replace(".sdf","")+"_out.sdf\" --log=\""+outpath+prot_file.replace(".pdb","")+"_"+lig_file.replace(".sdf","")+"_log.txt\" --seed=123456 --num_modes=10 --cpu=4")



























