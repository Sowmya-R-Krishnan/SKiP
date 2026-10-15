#Fix PDB files with missing atoms in backbone using PDBFixer and PyMOL

import os
import sys
import csv
import pandas as pd
import numpy as np
import pymol
from pymol import cmd

inpath = sys.argv[1]
outpath1 = sys.argv[2]
outpath2 = sys.argv[3]

#XXX: python 5_run_pdbfixer_and_pymol_v1.py "./Sample_output/FASPR_output/" "./Sample_output/Processed_PDBFixer/" "./Sample_output/Processed_PyMOL/"

for fname in os.listdir(inpath):
	if(fname.endswith(".pdb")):
		print(fname)
		
		os.system("pdbfixer "+inpath+fname+" --add-atoms=all --output "+outpath1+fname)
		
		cmd.load(outpath1+fname, "ref_chain")
		cmd.remove("hydrogens")
		cmd.save(outpath2+fname)
		cmd.reinitialize()
		
		#break




























