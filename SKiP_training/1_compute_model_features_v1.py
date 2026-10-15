#Program to compute enzyme and small molecule features for a given dataset

import os
import sys
sys.path.append("/home/sowmya/De_novo_DD/ChEMBL_dataset/Enzyme_engineering/BRENDA_database/Final_data/After_omissions/All_PDB_level_data/Structure_check/Model_datasets/")

import csv
import glob
import pandas as pd
import numpy as np
import Bio
import math
import subprocess
import shutil
import warnings
import rdkit
import msecif_v2
import pickle
from Bio import PDB, SeqUtils
from Bio.PDB import PDBIO, PDBParser, Select, NeighborSearch, SASA, ResidueDepth, NACCESS
from Bio.SeqUtils import seq1, IsoelectricPoint, ProtParam
from Bio.PDB.ResidueDepth import get_surface, residue_depth
from rdkit import Chem
from msecif_v2 import *
from collections import Counter
import pymol
from pymol import cmd

infile = sys.argv[1]
protpath = sys.argv[2]
docking_path = sys.argv[3]
outfile1 = sys.argv[4]  #Enzyme local features
outfile2 = sys.argv[5]  #Enzyme global features
outfile3 = sys.argv[6]  #Substrate features
tmp_path = sys.argv[7]
kegg_pathway_type = sys.argv[8]
naccess_binary = sys.argv[9]

#XXX: /usr/bin/python 1_compute_model_features_v1.py "./data/kcat_with_Km_WT_BS_data_logscale_final_v1.csv" "./Sample_structures/SKiD_kcat_archive/" "./Sample_structures/SKiD_kcat_archive/" "./Sample_output/Sample_local_feats_v1.csv" "./Sample_output/Sample_global_feats_v1.csv" "./Sample_output/Sample_substrate_feats_v1.csv" "./" "./data/Pathway_mapping_final.csv" <naccess binary path>

df = pd.read_csv(infile, sep="\t", header=0)
#print(df.columns)  #['EC_number', 'Substrate', 'UniProt_ID', 'Protein_file', 'Organism_name', 'kcat_value', 'Mutant', 'Mutation', 'pH', 'Temperature', 'References', 'Site_type', 'Substrate_SMILES', 'Mol_file', 'pkcat_value']

kegg_ptype = pd.read_csv(kegg_pathway_type, sep=",", header=0)
kegg_ptype.fillna(4, inplace=True)

#---------------------------------------------------------------FUNCTIONS------------------------------------------------------------------
#Function to compute substrate features using the mordred program
def compute_lig_feats(lig_file, tmp_path):	
	try:	
		os.system("python3 -m mordred -t sdf "+tmp_path+"tmp_out.sdf -o mordred.out -q -3")
		mol_df = pd.read_csv("mordred.out", sep=',', header=0)
		os.system("rm mordred.out")
		#print(mol_df)
		return mol_df
	except:
		print(lig_file+" failed.")
		return []
		
#Residue selection class
class ResidueSelect(Select):
	def __init__(self, accepted_resi, omit_resi):
		self.accepted_resi = accepted_resi
		self.omit_resi = omit_resi

	def accept_residue(self, resi):
		#print(atom.full_id)  #('6tf3', 0, 'A', ('H_3AT', 109, ' '), ('PB', ' '))
		if (resi in self.accepted_resi and resi not in self.omit_resi and not resi.id[0].startswith("H_")):
			return 1
		else:
			return 0
			
#Extract binding site residues given a PDB structure and ligand ID
def extract_binding_site_v2(ligand_id, pdb_file, pdb_path, outpath):
	with warnings.catch_warnings():
		warnings.simplefilter('ignore')
		pdb_id = pdb_file[0:-4]  #To omit *.pdb and *.ent extensions
		pdb = PDBParser().get_structure(pdb_id, pdb_path+pdb_file)
		io = PDBIO()
		io.set_structure(pdb)

		for model in pdb:
			for chain in model:
				for resi in chain:
					if(resi.id[0] == "H_"+ligand_id):
						center_atoms = Bio.PDB.Selection.unfold_entities(resi, 'A')
						alist = Bio.PDB.Selection.unfold_entities(pdb, 'A')
						search_center = NeighborSearch(alist)
						selected_resi = {res for center_atom in center_atoms for res in search_center.search(center_atom.coord, 6.5, 'R')}
						io.set_structure(pdb)
						io.save(outpath+"tmp_pocket.pdb", ResidueSelect(selected_resi, [resi]))
						
	selected_resi_final = [x for x in selected_resi if not x.id[0].startswith("H_")]
						
	return selected_resi_final
						
#Run NACCESS using the biopython interface instead of using shell commands
def get_asa_results(pdb_file, naccess_path):
	with warnings.catch_warnings():
		warnings.simplefilter('ignore')
		pdb = PDBParser().get_structure(pdb_file, pdb_file)
		io = PDBIO()
		io.set_structure(pdb)

		for model in pdb:
			output = NACCESS.NACCESS(model, pdb_file=pdb_file, naccess_binary=naccess_path, tmp_directory='/tmp')
			break  #Only for first model

		return output				
		
		
#Function to extract the PATTY atom type mapping for a list of PDB atoms
def extract_patty_types(pocket_file, patty_types):
	pdb_atoms = []
	with open(pocket_file) as f:
		for line in f.readlines():
			line = line.strip()
			pdb_atoms.append(line)
		
	patty_list = []

	for atom in pdb_atoms:
		resiname = atom[17:20].strip()
		atom_name = atom[12:16].strip()
		#To handle hydrogens present in crystal structures
		try:
			patty_type = patty_types[resiname][atom_name]
			patty_list.append(patty_type)
		except:
			continue

	return patty_list
	
#Function to translate a molecule to its centre of mass so that the principal moments are comparable
def translate_to_cofm(masses, xyz):
	# Position of centre of mass in original coordinates
	cofm = sum(masses[:,np.newaxis] * xyz) / np.sum(masses)
	# Transform to CofM coordinates and return
	xyz -= cofm
	return xyz

#Function to get the inertia matrix from masses and coordinates
def get_inertia_matrix(masses, xyz):
	# Moment of intertia tensor
	xyz = translate_to_cofm(masses, xyz)
	x, y, z = xyz.T
	Ixx = np.sum(masses * (y**2 + z**2))
	Iyy = np.sum(masses * (x**2 + z**2))
	Izz = np.sum(masses * (x**2 + y**2))
	Ixy = -np.sum(masses * x * y)
	Iyz = -np.sum(masses * y * z)
	Ixz = -np.sum(masses * x * z)
	I = np.array([[Ixx, Ixy, Ixz], [Ixy, Iyy, Iyz], [Ixz, Iyz, Izz]])
	return I

#Function to get the principal moments of inertia from the inertia matrix through diagonalization
def get_principal_moi(I):
	Ip = np.linalg.eigvals(I)
	# Sort and convert principal moments of inertia to SI (kg.m2)
	Ip.sort()
	return Ip

#Function to calculate euclidean distance between two points in Cartesian space of 3 dimensions
def euclidean_dist(x1, y1, z1, x2, y2, z2):
	return(math.sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2))
#------------------------------------------------------------------------------------------------------------------------------------------


out1 = open(outfile1, "w")
out2 = open(outfile2, "w")
out3 = open(outfile3, "w")

for i, row in df.iterrows():
	try:
		ph_val = float(row["pH"])
	except:
		ph_val = 7.0

	prot_file = row["Protein_file"].replace(".pdb","").replace("_mut","")+"_"+str(ph_val)
	lig_file = row["Mol_file"].replace(".sdf","")
	
	print(prot_file, lig_file)
	
	#-----------------------------------------------------------LIGAND FEATURES-------------------------------------------------------
	#XXX: As per the benchmarking of GNINA, the first docking pose will be used for feature calculations
	try:
		mol = Chem.SDMolSupplier(docking_path+row["Entry_ID"]+"/"+prot_file+"_"+lig_file+"_out.sdf")
	#Works for cases where the pH-dependent protonation has failed for some weird reasons
	except:
		print("Error: "+prot_file+"_"+lig_file+" is empty/missing")
		continue  #If ligand or protein feature calculations failed, the entire row is omitted
	
	for m in mol:
		conf = m.GetConformer()
		break  #Breaks the loop after first pose is saved
	
	writer = Chem.SDWriter(tmp_path+"tmp_out.sdf")	
	writer.write(m, confId=conf.GetId())
	writer.close()  #XXX: Very very important to close the file handle before passing on to mordred
	
	#Save the complex file also for local feature calculations
	cmd.load(protpath+row["Entry_ID"]+"/"+prot_file+".pdb", "mol1")
	cmd.load(tmp_path+"tmp_out.sdf", "mol2")
	cmd.save(tmp_path+"tmp_complex.pdb")
	cmd.reinitialize()
	
	mol_df = compute_lig_feats("tmp_out.sdf", tmp_path)
	mol_df["name"] = lig_file
	mol_df.fillna(0, inplace=True)

	if(i==0):
		cols = list(mol_df.columns)
		print("\t".join(cols), file=out3)
	
	out_str = []
	for rname, val in mol_df.iloc[0].items():
		out_str.append(str(val))
	
	#print(out_str)
	print("\t".join(out_str), file=out3)
	
	#--------------------------------------------------------------------------------------------------------------------------------

	#----------------------------------------------------GLOBAL ENZYME FEATURES------------------------------------------------------
	feats = {}
	feats["Protein_file"] = prot_file
	with warnings.catch_warnings():
		warnings.simplefilter('ignore')
		pdb = PDBParser().get_structure(prot_file, protpath+row["Entry_ID"]+"/"+prot_file+".pdb")
		io = PDBIO()
		io.set_structure(pdb)
		
		length = 0
		hetero_present = 0
		for model in pdb:
			for chain in model:
				seq = seq1(''.join(residue.resname for residue in chain if not residue.id[0].startswith("H_")))
				for resi in chain:
					if(not resi.id[0].startswith("H_")):  #Heteroatoms not counted as residues (might miss non-canonical residues?)
						length = length + 1
					else:
						hetero_present = hetero_present + 1
		
		feats["Length"] = length
		
		try:
			iep = IsoelectricPoint.IsoelectricPoint(seq)
			feats["IEP"] = iep.pi()
		except:
			print("Error: IEP calculation is not possible with BioPython - "+prot_file)
			continue
		
		if(hetero_present > 0):
			feats["Cofactor"] = 1
		else:
			feats["Cofactor"] = 0
		
		if(row["Mutant"]=="no"):  #WT enzymes
			feats["Mut"] = 0
			feats["Mut_type"] = 0
		else:  #Binding site mutations
			feats["Mut"] = 1
			feats["Mut_type"] = 1
			
		feats["EC_class"] = int(str(row["EC_number"]).split(".")[0])
			
		prot_params = ProtParam.ProteinAnalysis(seq)
		#print(dir(prot_params))  #['__class__', '__delattr__', '__dict__', '__dir__', '__doc__', '__eq__', '__format__', '__ge__', '__getattribute__', '__gt__', '__hash__', '__init__', '__init_subclass__', '__le__', '__lt__', '__module__', '__ne__', '__new__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__', '__str__', '__subclasshook__', '__weakref__', '_weight_list', 'amino_acids_content', 'amino_acids_percent', 'aromaticity', 'charge_at_pH', 'count_amino_acids', 'flexibility', 'get_amino_acids_percent', 'gravy', 'instability_index', 'isoelectric_point', 'length', 'molar_extinction_coefficient', 'molecular_weight', 'monoisotopic', 'protein_scale', 'secondary_structure_fraction', 'sequence']
		
		feats["Aromaticity"] = prot_params.aromaticity()
		if(row["pH"]!="-----"):
			feats["Charge"] = prot_params.charge_at_pH(float(row["pH"]))
		else:
			feats["Charge"] = prot_params.charge_at_pH(7.4)
			
		percent = prot_params.get_amino_acids_percent()
		for key, val in percent.items():
			feats["Percent_"+key] = val
		
		try:  #XXX: To handle: KeyError: 'U' - in biopython itself - Ex: PDB ID 6HN3
			feats["GRAVY_hydropathy"] = prot_params.gravy()
			feats["Instability_index"] = prot_params.instability_index()
			feats["Molar_extinction_coefficient_reduced"] = prot_params.molar_extinction_coefficient()[0]
			feats["Molar_extinction_coefficient_oxidized"] = prot_params.molar_extinction_coefficient()[1]
			feats["Molecular_weight"] = prot_params.molecular_weight()
			feats["Helix_fraction"] = prot_params.secondary_structure_fraction()[0]
			feats["Turn_fraction"] = prot_params.secondary_structure_fraction()[1]
			feats["Sheet_fraction"] = prot_params.secondary_structure_fraction()[2]
		except:
			print("Error: residue keyerror in BioPython - "+prot_file)
			feats["GRAVY_hydropathy"] = 0.0
			feats["Instability_index"] = 0.0
			feats["Molar_extinction_coefficient_reduced"] = 0.0
			feats["Molar_extinction_coefficient_oxidized"] = 0.0
			feats["Molecular_weight"] = 0.0
			feats["Helix_fraction"] = 0.0
			feats["Turn_fraction"] = 0.0
			feats["Sheet_fraction"] = 0.0
		
		#1 - Primary - Carbohydrate & Energy Metabolism; "2 - Primary - amino acids, fatty acids and nucleotides"; 3 - intermediate; 4 - secondary; 5 - x
		kegg_df = kegg_ptype[kegg_ptype["EC number"]==row["EC_number"]]
		try:
			feats["KEGG_pathway_type"] = int(kegg_df.iloc[0]["Pathway type"])
		except:
			feats["KEGG_pathway_type"] = 4 #Default is 4 - Secondary
		
		#feats["pKm_value"] = row["pKm_value"]  #TODO: To be uncommented only for kcat features
		
		if(i==0):
			cols_list = list(feats.keys())
			print("\t".join(cols_list), file=out2)
		
		vals = list(feats.values())
		vals_final = [str(x) for x in vals]
		print("\t".join(vals_final), file=out2)
	
	#----------------------------------------------------------------------------------------------------------------------
	
	#-----------------------------------------------------LOCAL ENZYME FEATURES--------------------------------------------
	l_feats = {}
	
	#Code adapted from https://github.com/koji11235/MSECIFv2/tree/main
	shells = [2.5, 4.5, 6.5, 8.5, 10.0]
	try:
		feature = msecif_v2.get_multi_shelled_ecif(protpath+row["Entry_ID"]+"/"+prot_file+".pdb", tmp_path+"tmp_out.sdf", shells)
	except:
		print("Error: PDB file is empty! - "+prot_file)
		continue
	
	columns = msecif_v2.make_feature_names_multi_shelled_ecif(shells)
	cols_final = ["ES_complex"]
	cols_final.extend(columns)
	feats_final = [prot_file+"_"+lig_file+"_out.sdf"]
	feats_final.extend(feature)
	for k,v in zip(cols_final, feats_final):
		l_feats[k] = v
	
	#Extract the binding site residues (without ligand or cofactor atoms) for local feature calculations
	bs_resi = extract_binding_site_v2("UNK", "tmp_complex.pdb", tmp_path, tmp_path)
	#print(bs_resi)
	#print(pfam_dict[prot_file[0:4]])
	
	if(len(bs_resi)>0):
		bs_resinum = []
		for resi in bs_resi:
			bs_resinum.append(resi.id[1])
		#print(bs_resinum)
		
		#Surface area features
		abs_sa = 0.0
		rel_sa = 0.0
		side_chain_abs = 0.0
		side_chain_rel = 0.0
		main_chain_abs = 0.0
		main_chain_rel = 0.0
		non_polar_abs = 0.0
		non_polar_rel = 0.0
		polar_abs = 0.0
		polar_rel = 0.0
		
		output = get_asa_results(protpath+row["Entry_ID"]+"/"+prot_file+".pdb", naccess_binary)
		#{('A', (' ', 1, ' ')): {'res_name': 'SER', 'all_atoms_abs': 127.59, 'all_atoms_rel': 109.5, 'side_chain_abs': 73.99, 'side_chain_rel': 94.7, 'main_chain_abs': 53.59, 'main_chain_rel': 139.6, 'non_polar_abs': 56.24, 'non_polar_rel': 115.8, 'all_polar_abs': 71.35, 'all_polar_rel': 105.0}
		for resi, asa in output.property_dict.items():
			if(resi[1][1] in bs_resinum):
				abs_sa = abs_sa + float(asa["all_atoms_abs"])
				rel_sa = rel_sa + float(asa["all_atoms_rel"])
				side_chain_abs = side_chain_abs + float(asa["side_chain_abs"])
				side_chain_rel = side_chain_rel + float(asa["side_chain_rel"])
				main_chain_abs = main_chain_abs + float(asa["main_chain_abs"])
				main_chain_rel = main_chain_rel + float(asa["main_chain_rel"])
				non_polar_abs = non_polar_abs + float(asa["non_polar_abs"])
				non_polar_rel = non_polar_rel + float(asa["non_polar_rel"])
				polar_abs = polar_abs + float(asa["all_polar_abs"])
				polar_rel = polar_rel + float(asa["all_polar_rel"])
				
		l_feats['abs_sa'] = abs_sa
		l_feats['rel_sa'] = rel_sa
		l_feats["side_chain_abs"] = side_chain_abs
		l_feats["side_chain_rel"] = side_chain_rel
		l_feats["main_chain_abs"] = main_chain_abs
		l_feats["main_chain_rel"] = main_chain_rel
		l_feats["non_polar_abs"] = non_polar_abs
		l_feats["non_polar_rel"] = non_polar_rel
		l_feats["polar_abs"] = polar_abs
		l_feats["polar_rel"] = polar_rel
		
		#Pharmacophore features: HBA (ACC), HBD (DON), Polar (POL), Anion (ANI), Cation (CAT), Hydrophobic (HYD), Others (OTH)
		#XXX: Based on PATTY atom typing in RDKit
		patty_types = {"GLY": {'N': 'CAT', 'CA': 'HYD', 'C': 'HYD', 'O': 'ACC'}, "ALA": {'N': 'CAT', 'CA': 'HYD', 'C': 'HYD', 'O': 'ACC', 'CB': 'HYD'}, "VAL": {'N': 'CAT', 'CA': 'HYD', 'C': 'HYD', 'O': 'ACC', 'CB': 'HYD', 'CG1': 'HYD', 'CG2': 'HYD'}, "LEU": {'N': 'CAT', 'CA': 'HYD', 'C': 'HYD', 'O': 'ACC', 'CB': 'HYD', 'CG': 'HYD', 'CD1': 'HYD', 'CD2': 'HYD'}, "ILE": {'N': 'CAT', 'CA': 'HYD', 'C': 'HYD', 'O': 'ACC', 'CB': 'HYD', 'CG1': 'HYD', 'CG2': 'HYD', 'CD1': 'HYD'}, "PHE": {'N': 'CAT', 'CA': 'HYD', 'C': 'HYD', 'O': 'ACC', 'CB': 'HYD', 'CG': 'HYD', 'CD1': 'HYD', 'CD2': 'HYD', 'CE1': 'HYD', 'CE2': 'HYD', 'CZ': 'HYD'}, "TYR": {'N': 'CAT', 'CA': 'HYD', 'C': 'HYD', 'O': 'ACC', 'CB': 'HYD', 'CG': 'HYD', 'CD1': 'HYD', 'CD2': 'HYD', 'CE1': 'HYD', 'CE2': 'HYD', 'CZ': 'OTH', 'OH': 'POL'}, "TRP": {'N': 'CAT', 'CA': 'HYD', 'C': 'HYD', 'O': 'ACC', 'CB': 'HYD', 'CG': 'HYD', 'CD1': 'HYD', 'CD2': 'HYD', 'CE2': 'OTH', 'CE3': 'HYD', 'NE1': 'DON', 'CZ2': 'HYD', 'CZ3': 'HYD', 'CH2': 'HYD'}, "MET": {'N': 'CAT', 'CA': 'HYD', 'C': 'HYD', 'O': 'ACC', 'CB': 'HYD', 'CG': 'HYD', 'SD': 'HYD', 'CE': 'HYD'}, "CYS": {'N': 'CAT', 'CA': 'HYD', 'C': 'HYD', 'O': 'ACC', 'CB': 'HYD', 'SG': 'HYD'}, "PRO": {'N': 'CAT', 'CA': 'HYD', 'C': 'HYD', 'O': 'ACC', 'CB': 'HYD', 'CG': 'HYD', 'CD': 'HYD'}, "SER": {'N': 'CAT', 'CA': 'HYD', 'C': 'HYD', 'O': 'ACC', 'CB': 'HYD', 'OG': 'POL'}, "THR": {'N': 'CAT', 'CA': 'HYD', 'C': 'HYD', 'O': 'ACC', 'CB': 'HYD', 'CG2': 'HYD', 'OG1': 'POL'}, "ASN": {'N': 'CAT', 'CA': 'HYD', 'C': 'HYD', 'O': 'ACC', 'CB': 'HYD', 'CG': 'OTH', 'ND2': 'DON', 'OD1': 'ACC'}, "GLN": {'N': 'CAT', 'CA': 'HYD', 'C': 'HYD', 'O': 'ACC', 'CB': 'HYD', 'CG': 'HYD', 'CD': 'OTH', 'NE2': 'DON', 'OE1': 'ACC'}, "ASP": {'N': 'CAT', 'CA': 'HYD', 'C': 'HYD', 'O': 'ACC', 'CB': 'HYD', 'CG': 'OTH', 'OD1': 'ANI', 'OD2': 'ANI'}, "GLU": {'N': 'CAT', 'CA': 'HYD', 'C': 'HYD', 'O': 'ACC', 'CB': 'HYD', 'CG': 'HYD', 'CD': 'OTH', 'OE1': 'ANI', 'OE2': 'ANI'}, "HIS": {'N': 'CAT', 'CA': 'HYD', 'C': 'HYD', 'O': 'ACC', 'CB': 'HYD', 'CG': 'OTH', 'CD2': 'HYD', 'ND1': 'POL', 'CE1': 'HYD', 'NE2': 'POL'}, "LYS": {'N': 'CAT', 'CA': 'HYD', 'C': 'HYD', 'O': 'ACC', 'CB': 'HYD', 'CG': 'HYD', 'CD': 'HYD', 'CE': 'HYD', 'NZ': 'CAT'}, "ARG": {'N': 'CAT', 'CA': 'HYD', 'C': 'HYD', 'O': 'ACC', 'CB': 'HYD', 'CG': 'HYD', 'CD': 'HYD', 'NE': 'CAT', 'CZ': 'OTH', 'NH1': 'CAT', 'NH2': 'CAT'}}

		#Pharmacophore mapping
		patty_list = extract_patty_types(tmp_path+"tmp_pocket.pdb", patty_types)
		patty_counter = Counter(patty_list)
		l_feats["HBA"] = patty_counter["ACC"]
		l_feats["HBD"] = patty_counter["DON"]
		l_feats["polar"] = patty_counter["POL"]
		l_feats["anion"] = patty_counter["ANI"]
		l_feats["cation"] = patty_counter["CAT"]
		l_feats["hydrophobic"] = patty_counter["HYD"]
		l_feats["others"] = patty_counter["OTH"]
		l_feats["total_heavy_atoms"] = len(patty_list)
		
		#Residue depth
		resi_depths = []
		with warnings.catch_warnings():
			warnings.simplefilter('ignore')
			pdb = PDBParser().get_structure(prot_file, "tmp_complex.pdb")
			io = PDBIO()
			io.set_structure(pdb)
			
			surface = get_surface(pdb[0], MSMS="/home/sowmya/PhD/RNA_drug_design_data/RNA_disease_information/Scripts/msms_i86_64Linux2_2.6.1/msms.x86_64Linux2.2.6.1")
			
			for model in pdb:
				for chain in model:
					for resi in chain:
						if(resi in bs_resi and not resi.id[0].startswith("H_")):
							resi_depths.append(residue_depth(chain[resi.id[1]], surface))
							
		l_feats["pocket_depth"] = np.mean(resi_depths)
		#print(np.mean(resi_depths))  #1U3W: 5.510
		
		#PMI descriptors
		mass_dict = {"H" :[1  , 1.00794], "He":[2  , 4.002602], "Li":[3  , 6.941], "Be":[4  , 9.012182], "B" :[5  , 10.811], "C" :[6  , 12.0107], "N" :[7  , 14.0067], "O" :[8  , 15.9994], "F" :[9  , 18.9984032], "Ne":[10 , 20.1797], "Na":[11 , 22.98976928], "Mg":[12 , 24.3050], "Al":[13 , 26.9815386], "Si":[14 , 28.0855], "P" :[15 , 30.973762], "S" :[16 , 32.065], "Cl":[17 ,  35.453], "Ar":[18 , 39.948], "K" :[19 , 39.0983], "Ca":[20 , 40.078], "Sc":[21 , 44.955912], "Ti":[22 , 47.867], "V" :[23 , 50.9415], "Cr":[24 , 51.9961], "Mn":[25 , 54.93845], "Fe":[26 ,55.933195], "Co":[27 , 58.933195], "Ni":[28 , 58.6934], "Cu":[29 , 63.546], "Zn":[30 , 65.39], "Ga":[31 , 69.723], "Ge":[32 , 72.64], "As":[33 , 74.92160], "Se":[34 , 78.96], "Br":[35 , 79.904], "Kr":[36 , 83.798], "Rb":[37 , 85.4678], "Sr":[38 , 87.62], "Y" :[39 , 88.90585], "Zr":[40 , 91.224], "Nb":[41 , 92.90638], "Mo":[42 , 95.94], "Tc":[43 , 97.9072], "Ru":[44 , 101.07], "Rh":[45 , 102.90550], "Pd":[46 , 106.42], "Ag":[47 , 107.8682], "Cd":[48 , 112.411], "In":[49 , 114.818], "Sn":[50 , 118.710], "Sb":[51 , 121.760], "Te":[52 , 127.60], "I" :[53 , 126.90447], "Xe":[54 , 131.293], "Cs":[55 , 132.9054519], "Ba":[56 , 137.327], "La":[57 , 138.90547], "Ce":[58 , 140.116], "Pr":[59 , 140.90765], "Nd":[60 , 144.242], "Pm":[61 , 145], "Sm":[62 , 150.36], "Eu":[63 , 151.964], "Gd":[64 , 157.25], "Tb":[65 , 158.92535], "Dy":[66 , 162.500], "Ho":[67 , 164.93032], "Er":[68 , 167.259], "Tm":[69 , 168.93421], "Yb":[70 , 173.04], "Lu":[71 , 174.967], "Hf":[72 , 178.49], "Ta":[73 , 180.94788], "W" :[74 , 183.84], "Re":[75 , 186.207], "Os":[76 , 190.23], "Ir":[77 , 192.217], "Pt":[78 , 195.084], "Au":[79 , 196.966569], "Hg":[80 , 200.59], "Tl":[81 , 204.3833], "Pb":[82 , 207.2], "Bi":[83 , 208.98040], "Po":[84 , 208.9824], "At":[85 , 209.9871], "Rn":[86 , 222.0176], "Fr":[87 , 223], "Ra":[88 , 226], "Ac":[89 , 227], "Th":[90 , 232.03806], "Pa":[91 , 231.03588], "U" :[92 , 238.02891], "Np":[93 , 238.8486], "Pu":[94 , 242.8798], "Am":[95 , 244.8594], "Cm":[96 , 246.911], "Bk":[97 , 248.9266], "Cf":[98 , 252.9578], "Es":[99 , 253.9656], "Fm":[100, 259.0046], "Md":[101, 260.0124], "No":[102, 261.0202], "Lr":[103, 264.0436]}

		masses = []
		coords = []
		with open(tmp_path+"tmp_pocket.pdb") as f:
			for line in f.readlines():
				line = line.strip()
				if(line.startswith("ATOM")):
					try:
						atom_type = line[12:16].strip()[0]
						#print(atom_type)
						masses.append(mass_dict[atom_type][1])
						coords.append([float(line[30:38].strip()), float(line[38:46].strip()), float(line[46:54].strip())])
					except:  #XXX: To handle errors in atom types (Ex: D in Phe164 of 5AE0)
						print("Error: Perdeuterated protein - "+prot_file)
						continue

		masses = np.asarray(masses)
		coords = np.asarray(coords)
		inertia_tensor = get_inertia_matrix(masses, coords)
		pmi = get_principal_moi(inertia_tensor)

		npr1 = np.round(pmi[0]/pmi[2], 5)
		npr2 = np.round(pmi[1]/pmi[2], 5)

		asphericity = np.round(0.5 * ((pmi[2]-pmi[0])**2 + (pmi[2]-pmi[1])**2 + (pmi[1]-pmi[0])**2)/(pmi[0]**2 + pmi[1]**2 + pmi[2]**2), 5)
		eccentricity = np.round(math.sqrt(pmi[2]**2 - pmi[0]**2)/pmi[2]**2, 5)
		spherocity_index = np.round((3*pmi[0])/(pmi[0]+pmi[1]+pmi[2]), 5)
		#rad_gyr = np.round(math.sqrt(2*3.14*((pmi[0]*pmi[1]*pmi[2])/3)/(np.sum(masses))), 5)
		inertial_shape_factor = np.round(pmi[1]/(pmi[0]*pmi[2]), 5)

		#print(npr1, npr2, asphericity, eccentricity, spherocity_index, inertial_shape_factor)
		l_feats["pmi_1"] = pmi[0]
		l_feats["pmi_2"] = pmi[1]
		l_feats["pmi_3"] = pmi[2]
		l_feats["npr1"] = npr1
		l_feats["npr2"] = npr2
		l_feats["asphericity"] = asphericity
		l_feats["eccentricity"] = eccentricity
		l_feats["spherocity_index"] = spherocity_index
		l_feats["inertial_shape_factor"] = inertial_shape_factor
		
		os.system("rm tmp_complex.pdb tmp_out.sdf tmp_pocket.pdb")
		
		if(i==0):
			cols_list = list(l_feats.keys())
			print("\t".join(cols_list), file=out1)
		
		vals = list(l_feats.values())
		vals_final = [str(x) for x in vals]
		print("\t".join(vals_final), file=out1)			
	else:
		continue		
	
	#----------------------------------------------------------------------------------------------------------------------------
	#if(i>9):
		#break

out1.close()
out2.close()
out3.close()


























