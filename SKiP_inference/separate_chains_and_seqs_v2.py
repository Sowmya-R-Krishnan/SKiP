#Program to save each chain of a PDB structure and its fasta sequence (without missing residues) using BioPython

import sys
import csv
import pandas as pd
import numpy as np
import Bio
import warnings

from Bio import PDB, SeqUtils
from Bio.PDB import PDBIO, PDBParser, Select
from Bio.SeqUtils import seq1
from collections import defaultdict

input = sys.argv[1]  #PDB ID
inpath = sys.argv[2]  #Path to the full PDB structure
mutstr = sys.argv[3]  #Mutations to be introduced in the PDB-derived sequence - Ex: A:H110R;A:P23G
outpath = sys.argv[4]  #Path to write the chain-wise PDB and Fasta files

#print(input, inpath, mutstr, outpath)

#pdb_id = input.split("_")[0]
pdb_id=input
mutations = mutstr.split("+")  #XXX: Use single-letter residue format instead of 3-letter

#Residue selection class
class ResidueSelect(Select):
# # A sample method  def fun(self): where fun describes some steps self is part of default syntax
	def __init__(self, accepted_resi, omit_resi):
		self.accepted_resi = accepted_resi
		self.omit_resi = omit_resi

	def accept_residue(self, resi):
		#print(atom.full_id)  #('6tf3', 0, 'A', ('H_3AT', 109, ' '), ('PB', ' '))
		if (resi in self.accepted_resi and resi not in self.omit_resi):
			return 1
		else:
			return 0
#https://www.tutorialspoint.com/what-is-the-use-of-the-with-statement-in-python
with warnings.catch_warnings():
	warnings.simplefilter('ignore')
	#print(pdb_id)
	#All warnings are ignored by setting the first argument action of warnings.simplefilter() to 'ignore'
	pdb = PDBParser().get_structure(pdb_id, inpath+pdb_id+".pdb") #PDBParser module to get structure
	io = PDBIO()
	io.set_structure(pdb)#The input structure
	
	chains = {chain.id:seq1(''.join(residue.resname for residue in chain if not residue.id[0].startswith("H_"))) for chain in pdb.get_chains()}#All the chains in the pdb can be retrieved with get_chains
	#from Bio.SeqUtils import seq1
#one-liner to extract the amino acid sequence of a specific chain from a PDB file
	#print(chains)
	chain_residict = defaultdict(dict)
	for chain in pdb.get_chains():
		chain_residict[chain.id] = {}
		#print(chain_residict[chain.id])
		for residue in chain:
			if(not residue.id[0].startswith("H_") and not residue.id[0].startswith("W")):
				residue_id = str(residue.id[1])
				if(residue.id[2]!=" "):
					residue_id = residue_id+str(residue.id[2])  #To handle residue insertions like "355A" - Ex: 3QYA
				chain_residict[chain.id][residue_id] = seq1(residue.resname)
			#print(chain_residict[chain.id][residue.id[1]])
	
	for model in pdb:
		for chain in model:
			selected_resi = []
			for resi in chain:
				selected_resi.append(resi)
					
			io.set_structure(pdb)
			io.save(outpath+pdb_id+"_"+str(chain.id)+".pdb", ResidueSelect(selected_resi, []))
			
			#WT sequence FASTA file
			with open(outpath+pdb_id+"_"+str(chain.id)+".fasta", "w") as out:
				#TODO 1: Remove X corresponding to water molecules, modified residues, HETATM records
				seq = chains[chain.id].replace("X","")
				#seq = str.lower(seq)
				#print (seq)
				#print(">"+pdb_id+"_"+chain.id, file=out)
				print(seq, file=out)
				
			#Mutant sequence FASTA file
			mut_seq = ""
			counter = 1
			for mut in mutations:
				chain_id = mut.split(":")[0]
				if(chain_id==chain.id):
					chain_residict_copy = chain_residict[chain_id]
					#print(chain_residict[chain_id])
					if(counter == 1):
						for k in chain_residict_copy.keys():
							chain_residict_copy[k]=str.lower(chain_residict_copy[k])
						#print(chain_residict_copy[k])
					res_id = int(mut.split(":")[1][1:-1])
					#print (res_id)
					mutant_resi = mut.split(":")[1][-1]
					#print (mutant_resi)
					chain_residict_copy[str(res_id)] = mutant_resi
					mut_seq = "".join(list(chain_residict_copy.values())).replace("x","")
					counter = 2
				else:
					continue
			
			if(mut_seq!=""):	
				with open(outpath+pdb_id+"_"+str(chain.id)+"_"+str(mutstr)+"_mut.fasta", "w") as out:
					print(mut_seq, file=out)
				
			out.close()
	
print("Done")

























