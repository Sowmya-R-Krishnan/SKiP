#Program to extract mutation datasets for running FASPR automation code

import sys
import csv
import pandas as pd
import numpy as np
from collections import defaultdict

infile = sys.argv[1]
outfile = sys.argv[2]

#XXX: python 3_prepare_FASPR_input_v1.py "./Sample_output/kcat_test_set_PDB_muts_v1.csv" "./Sample_output/kcat_FASPR_input_v1.csv"

mut_mapped_df = pd.read_csv(infile, sep="\t", header=0)

three_to_one = {"ALA":"A", "GLY":"G", "VAL":"V", "LEU":"L", "ILE":"I", "PHE":"F", "TRP":"W", "TYR":"Y", "HIS":"H", "ARG":"R", "LYS":"K", "ASN":"N", "GLN":"Q", "ASP":"D", "GLU":"E", "SER":"S", "CYS":"C", "THR":"T", "MET":"M", "PRO":"P", "OCS":"C", "CSO":"C", "MSE":"M", "MHS":"H", "CME":"C", "CSD":"C", "TSY":"C", "CSS":"C", "SCY":"C", "PTR":"Y", "DDZ":"A", "JJJ":"C", "MHO":"M", "TPO":"T", "SEP":"S", "KPI":"K", "SLZ":"K", "KPF":"K", "KPY":"K", "74P":"K", "VPV":"K", "KGC":"K", "KYQ":"K", "LYF":"K", "CIR":"R", "SNC":"C", "CSP":"C", "LLP":"K", "KCX":"K"}

out = open(outfile, "w")
print("affin.pdb_id\taffin.chain\tmutation", file=out)

for i, row in mut_mapped_df.iterrows():
	if(row["Mutant"]=="yes"):
		pdb = row["PDB_ID"]
		brenda_mut = row["Mutations"].split("/")
		PDB_mut = row["PDB_muts"].split("/")
		
		mut_str = []
		for bmut, pmut in zip(brenda_mut, PDB_mut):  #TODO: Assumes that both mutation lists follow the same order! Make sure its valid!!
			pdb_mut_new = three_to_one[pmut.split("_")[0]]+pmut.split("_")[1]+bmut[-1]
			mut_str.append(pdb_mut_new)
		
		print(str(pdb)+"\t"+str(row["Chain"])+"\t"+"/".join(mut_str), file=out)
	
out.close()


























