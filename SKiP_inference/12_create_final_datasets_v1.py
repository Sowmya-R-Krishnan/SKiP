#Program to separate the features for kcat and Km model checks

import sys
import csv
import pandas as pd
import numpy as np

ref_file = sys.argv[1]  #File with the list of features selected for the model training as header
main_dataset = sys.argv[2]
local_feats = sys.argv[3]
global_feats = sys.argv[4]
sub_feats = sys.argv[5]
rdm_feats = sys.argv[6]
outfile = sys.argv[7]

#XXX: python 12_create_final_datasets_v1.py "../SKiP_training/Sample_output/Random_split/kcat_X_test.csv" "./data/kcat_independent_test_set_v1.csv" "./Sample_output/Features/kcat_local_feats_v1.csv" "./Sample_output/Features/kcat_global_feats_v1.csv" "./Sample_output/Features/kcat_substrate_feats_v1.csv" "./Sample_output/Features/kcat_RDM_fingerprints_v1.csv" "./Sample_output/Features/kcat_independent_test_set_feats_final_v1.csv"

ref_file = pd.read_csv(ref_file, sep="\t", header=0)
nomit = 3
feats_final = list(ref_file.columns)

main_df = pd.read_csv(main_dataset, sep="\t", header=0)
local_df = pd.read_csv(local_feats, sep="\t", header=0)
global_df = pd.read_csv(global_feats, sep="\t", header=0)
sub_df = pd.read_csv(sub_feats, sep="\t", header=0)
rdm_df = pd.read_csv(rdm_feats, sep="\t", header=0)

col_list = []
col_list.extend(list(global_df.columns))
col_list.insert(1, list(local_df.columns)[0])
col_list.extend(list(local_df.columns)[1:])
col_list.insert(2, list(sub_df.columns)[0])
col_list.extend(list(sub_df.columns)[1:])
col_list.extend(list(rdm_df.columns[2:]))
merged_df = pd.DataFrame(columns=col_list)

three_to_one = {"ALA":"A", "GLY":"G", "VAL":"V", "LEU":"L", "ILE":"I", "PHE":"F", "TRP":"W", "TYR":"Y", "HIS":"H", "ARG":"R", "LYS":"K", "ASN":"N", "GLN":"Q", "ASP":"D", "GLU":"E", "SER":"S", "CYS":"C", "THR":"T", "MET":"M", "PRO":"P", "OCS":"C", "CSO":"C", "MSE":"M", "MHS":"H", "CME":"C", "CSD":"C", "TSY":"C", "CSS":"C", "SCY":"C", "PTR":"Y", "DDZ":"A", "JJJ":"C", "MHO":"M", "TPO":"T", "SEP":"S", "KPI":"K", "SLZ":"K", "KPF":"K", "KPY":"K", "74P":"K", "VPV":"K", "KGC":"K", "KYQ":"K", "LYF":"K", "CIR":"R", "SNC":"C", "CSP":"C", "LLP":"K", "KCX":"K"}

for i, row in main_df.iterrows():
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
	
	prot = prot_file+"_"+str(ph_val)
	lig = row["Mol_ID"]
	
	g_sub_df = global_df[global_df["Protein_file"]==prot]
	lig_sub_df = sub_df[sub_df["name"]==lig]
	l_sub_df = local_df[local_df["ES_complex"]==prot+"_"+lig+"_out.sdf"]
	rdm_sub_df = rdm_df[rdm_df["ES_complex"]==prot+"_"+lig]
	
	if(len(list(g_sub_df.index))==0):
		g_sub_df = global_df[global_df["Protein_file"]==prot]
		l_sub_df = local_df[local_df["ES_complex"]==prot+"_"+lig+"_out.sdf"]
		
		if(len(g_sub_df.index)==0):
			print(prot, lig)
			row_df = pd.DataFrame([row])
			failed_df = pd.concat([failed_df, row_df], ignore_index=True)
			continue  #Move on to the next row
			
	elif(len(list(l_sub_df.index))==0 or len(list(lig_sub_df.index))==0):
		continue  #Move on to the next row
		
	row_dict = {}
	g_dict = pd.DataFrame([g_sub_df.iloc[0]]).to_dict(orient='list')
	l_dict = pd.DataFrame([l_sub_df.iloc[0]]).to_dict(orient='list')
	lig_dict = pd.DataFrame([lig_sub_df.iloc[0]]).to_dict(orient='list')
	rdm_dict = pd.DataFrame([rdm_sub_df.iloc[0]]).to_dict(orient='list')
	row_dict.update(g_dict)
	row_dict.update(l_dict)
	row_dict.update(lig_dict)
	row_dict.update(rdm_dict)
	#row_dict["pKm_value"] = [row["pKm_value"]]  #TODO: Uncomment only when computing features for the kcat dataset
	datapoint = pd.DataFrame.from_dict(row_dict)
	datapoint.fillna(0, inplace=True)
	
	merged_df = pd.concat([merged_df, datapoint], ignore_index=True)

sub_df = merged_df[feats_final]
sub_df.to_csv(outfile, sep="\t", header=True, index=False)















