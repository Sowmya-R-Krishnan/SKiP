#Program to extract reaction counts per EC based on the mapped RDM

import os
import sys
import csv
import pandas as pd
import numpy as np
from collections import defaultdict

infile = sys.argv[1]
rdm = sys.argv[2]
cof_rdm = sys.argv[3]
outfile1 = sys.argv[4]
outfile2 = sys.argv[5]

#XXX: python 10_get_reaction_frequency_v1.py "./Sample_output/kcat_MetaCyc_rxns_v1.csv" "../SKiP_training/data/uniqueRDM_db.csv" "../SKiP_training/data/cofactor_EC_v02 1.csv" "./Sample_output/kcat_RDMs_unmapped_v1.csv" "./Sample_output/kcat_RDMs_mapped_v1.csv"

df = pd.read_csv(infile, sep="\t", header=0)
rdm_df = pd.read_csv(rdm, sep="\t", header=0)
cof_rdm_df = pd.read_csv(cof_rdm, sep="\t", header=0)

ec_rdm_counts = defaultdict(list)
ec_rxn_rdm_map = defaultdict(list)
rdm_sub_df = rdm_df[rdm_df["Matched-part"]=="RDM"]

for i, row in rdm_sub_df.iterrows():
	ec_list = row["EC-wise-count"].split("||")
	rxn_list = row["EC-wise-reactionID"].split("||")
	for ec in ec_list:
		cof_check = cof_rdm_df[(cof_rdm_df["sorted_RDM"]==row["RDM"]) & (cof_rdm_df["Co-factor_list_ec"].str.contains(ec.split(":")[0].replace("EC-","")))]
		
		if(len(cof_check.index)==0):
			ec_rdm_counts[ec.split(":")[0]].append((row["RDM"], int(ec.split(":")[1])))
			
	for ec in rxn_list:
		ec_num = ec.split(":")[0]
		rxn_ids = ec.split(":")[-1].split(",")
		for rid in rxn_ids:
			if(row["RDM"] not in ec_rxn_rdm_map[(ec_num, rid)]):
				ec_rxn_rdm_map[(ec_num, rid)].append(row["RDM"])

sub_df = df[df["TC_value"]<0.5]  #Unmapped substrates with low similarity to MetaCyc reactants and products

out = open(outfile1, "w")
print("Entry_ID\tEC_number\tPDB_ID\tpH\tSubstrate_SMILES\tMol_ID\tRDM_1\tRDM_2\tRDM_3", file=out)

sub_ec_rdm_counts = defaultdict(list)
for i, row in sub_df.iterrows():
	ec_parts = row["EC_number"].split(".")
	new_ec = "EC-"+ec_parts[0]+"."+ec_parts[1]+"."+ec_parts[2]
	
	rdm_list = ec_rdm_counts[new_ec]
	sorted_rdm = sorted(rdm_list, key=lambda x: x[1], reverse=True)
	
	if(len(sorted_rdm)>=3):
		print(row["Entry_ID"]+"\t"+row["EC_number"]+"\t"+row["PDB_ID"]+"\t"+str(row["pH"])+"\t"+row["Substrate_SMILES"]+"\t"+str(row["Mol_ID"])+"\t"+str(sorted_rdm[0][0])+"\t"+str(sorted_rdm[1][0])+"\t"+str(sorted_rdm[2][0]), file=out)
	elif(len(sorted_rdm)==2):
		print(row["Entry_ID"]+"\t"+row["EC_number"]+"\t"+row["PDB_ID"]+"\t"+str(row["pH"])+"\t"+row["Substrate_SMILES"]+"\t"+str(row["Mol_ID"])+"\t"+str(sorted_rdm[0][0])+"\t"+str(sorted_rdm[1][0])+"\t-", file=out)
	elif(len(sorted_rdm)==1):
		print(row["Entry_ID"]+"\t"+row["EC_number"]+"\t"+row["PDB_ID"]+"\t"+str(row["pH"])+"\t"+row["Substrate_SMILES"]+"\t"+str(row["Mol_ID"])+"\t"+str(sorted_rdm[0][0])+"\t-\t-", file=out)
	else:
		print(row["Entry_ID"], row["EC_number"], row["Mol_ID"])
	
	#break

out.close()

sub_df = df[df["TC_value"]>=0.5]  #Mapped substrates with low similarity to MetaCyc reactants and products

out = open(outfile2, "w")
print("Entry_ID\tEC_number\tPDB_ID\tpH\tSubstrate_SMILES\tMol_ID\tMetaCyc_reaction_ID\tRDM_list", file=out)

sub_ec_rdm_counts = defaultdict(list)
for i, row in sub_df.iterrows():
	ec_parts = row["EC_number"].split(".")
	new_ec = "EC-"+ec_parts[0]+"."+ec_parts[1]+"."+ec_parts[2]
	rdm_list = ec_rxn_rdm_map[(new_ec, row["Reaction_ID"])]
	
	if(len(rdm_list)>0):
		print(row["Entry_ID"]+"\t"+row["EC_number"]+"\t"+row["PDB_ID"]+"\t"+str(row["pH"])+"\t"+row["Substrate_SMILES"]+"\t"+str(row["Mol_ID"])+"\t"+str(row["Reaction_ID"])+"\t"+str("||".join(rdm_list)), file=out)
	else:
		print(row["Entry_ID"], row["EC_number"], row["Mol_ID"], row["Reaction_ID"])
	
	#break

out.close()


















