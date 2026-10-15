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

#XXX: python 4_compute_R_part_fingerprint_v1.py "./data/kcat_with_Km_WT_BS_data_logscale_final_v1.csv" "./Sample_output/kcat_dataset_RDMs_mapped_sample.csv" "./Sample_output/kcat_dataset_RDMs_unmapped_sample.csv" "./data/Atom_types.xlsx" "./Sample_output/kcat_sample_RDM_fingerprints_v1.csv" "./data/In_house_RDMs.csv"

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
		if((row["EC_number"], row["Mol_file"]) in list(done.keys())):
			idx = done[(row["EC_number"], row["Mol_file"])]
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
			
		prot_file = row["Protein_file"].replace(".pdb","").replace("_mut","")		
		es_complex = prot_file+"_"+row["Mol_file"].replace(".sdf","_out.sdf")
		counts.insert(0, row["Entry_ID"])
		counts.insert(1, es_complex)
		row_dict = dict((col, [x]) for col, x in zip(cols_final, counts))
		row_df = pd.DataFrame.from_dict(row_dict)
		final_df = pd.concat([final_df, row_df], ignore_index=True)
		
		if(redundant==0):
			done[(row["EC_number"], row["Mol_file"])] = i
	else:
		print(row["Entry_ID"])

final_df.to_csv(outfile, sep="\t", header=True, index=False)




























