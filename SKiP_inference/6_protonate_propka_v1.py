#Protonate the protein based on pH from BRENDA using PDB2PQR and PROPKA tools

import os
import sys
import csv
import pandas as pd
import numpy as np
import pymol
from pymol import cmd

infile = sys.argv[1]
pdbpath = sys.argv[2]
outpath = sys.argv[3]
logpath = sys.argv[4]

#XXX: python 6_protonate_propka_v1.py "./data/kcat_independent_test_set_v1.csv" "./Sample_output/Processed_PyMOL/" "./Sample_output/Docking_ready/" "./Sample_output/Propka_logs/"

df = pd.read_csv(infile, sep="\t", header=0)
fnames = os.listdir(pdbpath)

three_to_one = {"ALA":"A", "GLY":"G", "VAL":"V", "LEU":"L", "ILE":"I", "PHE":"F", "TRP":"W", "TYR":"Y", "HIS":"H", "ARG":"R", "LYS":"K", "ASN":"N", "GLN":"Q", "ASP":"D", "GLU":"E", "SER":"S", "CYS":"C", "THR":"T", "MET":"M", "PRO":"P", "OCS":"C", "CSO":"C", "MSE":"M", "MHS":"H", "CME":"C", "CSD":"C", "TSY":"C", "CSS":"C", "SCY":"C", "PTR":"Y", "DDZ":"A", "JJJ":"C", "MHO":"M", "TPO":"T", "SEP":"S", "KPI":"K", "SLZ":"K", "KPF":"K", "KPY":"K", "74P":"K", "VPV":"K", "KGC":"K", "KYQ":"K", "LYF":"K", "CIR":"R", "SNC":"C", "CSP":"C", "LLP":"K", "KCX":"K"}

done = []
for i, row in df.iterrows():
	prot_file = row["PDB_ID"]
	if(row["Mutant"]=="yes"):
		brenda_mut = row["Mutations"].split("/")
		PDB_mut = row["PDB_muts"].split("/")
		
		mut_str = []
		for bmut, pmut in zip(brenda_mut, PDB_mut):  #TODO: Assumes that both mutation lists follow the same order! Make sure its valid!!
			pdb_mut_new = three_to_one[pmut.split("_")[0]]+pmut.split("_")[1]+bmut[-1]
			mut_str.append(row["Chain"]+":"+pdb_mut_new)
			
		prot_file = row["PDB_ID"]+"_"+row["Chain"]+"_"+"+".join(mut_str)+"_mut"
		
	out_file = prot_file.replace(".pdb","").replace("_mut","")  #Covers both mutant and wild-type structures
	pdb_id = prot_file[0:4]
	
	try:
		ph_val = float(row["pH"])
	except:
		ph_val = 7.0
	
	#Check if PDB file is present
	if((out_file+".pdb" in fnames or out_file+"_mut.pdb" in fnames) and not os.path.exists(outpath+out_file+"_"+str(ph_val)+".pdb")):
		#Remove ligand atoms and save the protein file
		het_records = []
		prot_records = []
		with open(pdbpath+prot_file+".pdb") as f:
			for line in f.readlines():
				line = line.strip()
				if(line.startswith("HETATM")):
					het_records.append(line)
				else:
					prot_records.append(line)
					
		temp_out = open(outpath+"tmp.pdb", "w")
		for line in prot_records:
			print(line, file=temp_out)
		temp_out.close()
		
		#Run PDB2PQR and PROPKA tools for protonating as per pH
		os.system("pdb2pqr "+outpath+"tmp.pdb "+outpath+"tmp.pqr --ff=CHARMM --pdb-output="+outpath+"tmp_1.pdb --with-ph="+str(ph_val)+" --titration-state-method=propka -o "+str(ph_val)+" -q")
		os.system("rm "+outpath+"tmp.pdb")
		os.system("rm "+outpath+"tmp.pqr")
		os.system("mv "+outpath+"tmp.log "+logpath+out_file+"_"+str(ph_val)+".log")
		
		#Append substrate and cofactor atoms back to the protonated PDB file
		mod_prot = open(outpath+out_file+"_"+str(ph_val)+".pdb", "w")
		with open(outpath+"tmp_1.pdb") as f:
			for line in f.readlines():
				line = line.strip()
				if(line!="END"):
					print(line, file=mod_prot)
				else:
					for l in het_records:
						print(l, file=mod_prot)
					print("END", file=mod_prot)
					
		mod_prot.close()
		os.system("rm "+outpath+"tmp_1.pdb")
		done.append(out_file)		
	else:
		continue
	
	print(out_file)

	#break

























