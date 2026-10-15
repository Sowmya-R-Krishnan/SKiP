#Program to compute RDM R-part fingerprint for all SKiD entries

import os
import sys
import csv
import pandas as pd
import numpy as np
from collections import Counter

dataset = sys.argv[1]
mapped = sys.argv[2]
unmapped = sys.argv[3]
atom_types = sys.argv[4]
outfile = sys.argv[5]
addn_data = sys.argv[6]

#XXX: python 11_compute_R_part_fingerprint_v1.py "./data/kcat_independent_test_set_v1.csv" "./Sample_output/kcat_RDMs_mapped_v1.csv" "./Sample_output/kcat_RDMs_unmapped_v1.csv" "../SKiP_training/data/Atom_types.xlsx" "./Sample_output/Features/kcat_RDM_fingerprints_v1.csv" "../SKiP_training/data/In_house_RDMs.csv" 

df = pd.read_csv(dataset, sep="\t", header=0)
mapped_df = pd.read_csv(mapped, sep="\t", header=0)
unmapped_df = pd.read_csv(unmapped, sep="\t", header=0)
kegg_atoms = pd.read_excel(atom_types, header=0)
addn_df = pd.read_csv(addn_data, sep=",", header=0)

bits = list(kegg_atoms["Atom type"])
bit_encoding = dict((bit, i+1) for i, bit in enumerate(bits))
print("No. of KEGG atom types in RDMs: "+str(len(bits)))

cols_final = ["Entry_ID", "ES_complex"]
cols_final.extend(bits)
final_df = pd.DataFrame(columns=cols_final)

three_to_one = {"ALA":"A", "GLY":"G", "VAL":"V", "LEU":"L", "ILE":"I", "PHE":"F", "TRP":"W", "TYR":"Y", "HIS":"H", "ARG":"R", "LYS":"K", "ASN":"N", "GLN":"Q", "ASP":"D", "GLU":"E", "SER":"S", "CYS":"C", "THR":"T", "MET":"M", "PRO":"P", "OCS":"C", "CSO":"C", "MSE":"M", "MHS":"H", "CME":"C", "CSD":"C", "TSY":"C", "CSS":"C", "SCY":"C", "PTR":"Y", "DDZ":"A", "JJJ":"C", "MHO":"M", "TPO":"T", "SEP":"S", "KPI":"K", "SLZ":"K", "KPF":"K", "KPY":"K", "74P":"K", "VPV":"K", "KGC":"K", "KYQ":"K", "LYF":"K", "CIR":"R", "SNC":"C", "CSP":"C", "LLP":"K", "KCX":"K"}

done = {}
for i, row in df.iterrows():
	if(row["Entry_ID"] in list(mapped_df["Entry_ID"])):
		sub_df = mapped_df[mapped_df["Entry_ID"]==row["Entry_ID"]]
	else:
		sub_df = unmapped_df[unmapped_df["Entry_ID"]==row["Entry_ID"]]
		
	if(len(sub_df.index)==0):
		sub_df = addn_df[(addn_df["EC_number"]==row["EC_number"]) & (addn_df["Substrate"]==row["Substrate"])]

	redundant = 0
	if(len(sub_df.index)==0):
		if((row["EC_number"], row["Mol_ID"]) in list(done.keys())):
			idx = done[(row["EC_number"], row["Mol_ID"])]
			entry_id = df.iloc[idx]["Entry_ID"]
			redundant = 1
	
			sub_df = mapped_df[mapped_df["Entry_ID"]==entry_id]
			if(len(sub_df.index)==0):
				sub_df = unmapped_df[unmapped_df["Entry_ID"]==entry_id]
	
	if(len(sub_df.index)>0):
		try:
			rdms = sub_df.iloc[0]["RDM_list"].split("||")
		except:
			rdms = []
			rdms.append(sub_df.iloc[0]["RDM_1"])
			rdms.append(sub_df.iloc[0]["RDM_2"])
			rdms.append(sub_df.iloc[0]["RDM_3"])
		
		r_parts = [r.split(":")[0] for r in rdms if r!="-"]
		counts = [0]*len(bits)
		for r_part in r_parts:
			atype1 = r_part.split("-")[0].strip()
			atype2 = r_part.split("-")[1].strip()
			counts[bit_encoding[atype1]-1] += 1
			counts[bit_encoding[atype2]-1] += 1
		
		try:
			ph_val = float(row["pH"])
		except:
			ph_val = 7.0
			
		prot_file = row["PDB_ID"]
		if(row["Mutant"]=="yes"):
			brenda_mut = row["Mutations"].split("/")
			PDB_mut = row["PDB_muts"].split("/")
			
			mut_str = []
			for bmut, pmut in zip(brenda_mut, PDB_mut):  #TODO: Assumes that both mutation lists follow the same order! Make sure its valid!!
				pdb_mut_new = three_to_one[pmut.split("_")[0]]+pmut.split("_")[1]+bmut[-1]
				mut_str.append(row["Chain"]+":"+pdb_mut_new)
				
			prot_file = row["PDB_ID"]+"_"+row["Chain"]+"_"+"+".join(mut_str)
		
		es_complex = prot_file+"_"+str(ph_val)+"_"+row["Mol_ID"].replace(".sdf","_out.sdf")
		counts.insert(0, row["Entry_ID"])
		counts.insert(1, es_complex)
		row_dict = dict((col, [x]) for col, x in zip(cols_final, counts))
		row_df = pd.DataFrame.from_dict(row_dict)
		final_df = pd.concat([final_df, row_df], ignore_index=True)
		
		if(redundant==0):
			done[(row["EC_number"], row["Mol_ID"])] = i
	else:
		print(row["Entry_ID"])

final_df.to_csv(outfile, sep="\t", header=True, index=False)




























