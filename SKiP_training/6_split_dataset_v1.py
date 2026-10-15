#Program to perform different random and stratified splits on the dataset

import os
import sys
import csv
import pandas as pd
import numpy as np
from sklearn import model_selection
from sklearn.model_selection import train_test_split
from collections import Counter, defaultdict

infile = sys.argv[1]
split_type = sys.argv[2]
act_type = sys.argv[3]
ec_map = sys.argv[4]
outpath = sys.argv[5]

#XXX: Random split = python 6_split_dataset_v1.py "./Sample_output/kcat_sample_final_feats_v1.csv" "random" "kcat" "./data/kcat_with_Km_WT_BS_data_logscale_final_v1.csv" "./Sample_output/Random_split/"

#XXX: EC-wise split = python 6_split_dataset_v1.py "./Sample_output/kcat_sample_final_feats_v1.csv" "strat_ec_1" "kcat" "./data/kcat_with_Km_WT_BS_data_logscale_final_v1.csv" "./Sample_output/EC_wise_split/"

'''
Train-val-test split types
==========================
a) random - Random 80-10-10 split
b) strat_ec_1 - Stratified split 1 on enzymes - data points of each EC number appears in 80:20 ratio in training and test datasets
'''

df = pd.read_csv(infile, sep="\t", header=0)
X_dataset = df[list(df.columns)[:-1]]
Y_dataset = df[list(df.columns)[-1]]

if(split_type not in ["random", "strat_ec_1"]):
	print('Split type not supported. Must be one of {"random", "strat_ec_1"}.')
else:
	if(split_type=="random"):
		random_state = 42  #XXX: Very very important for reproducibility of the split
		X_train, X_test, Y_train, Y_test = train_test_split(X_dataset, Y_dataset, test_size=0.1, random_state=random_state, shuffle=True)
		print(len(X_train), len(X_test), len(Y_train), len(Y_test))
		
		X_train.to_csv(outpath+act_type+"_X_train.csv", sep="\t", header=True, index=False)
		Y_train.to_csv(outpath+act_type+"_Y_train.csv", sep="\t", header=True, index=False)
		X_test.to_csv(outpath+act_type+"_X_test.csv", sep="\t", header=True, index=False)
		Y_test.to_csv(outpath+act_type+"_Y_test.csv", sep="\t", header=True, index=False)
		
		print("Random shuffled split done.")
		
	elif(split_type=="strat_ec_1"):
		ec_df = pd.read_csv(ec_map, sep="\t", header=0)  #Required since the features file does not have EC numbers mapped
		ec_dict = {}
		prot_dict = defaultdict(list)
		for i, row in ec_df.iterrows():
			ec_dict[row["Parent_PDB_structure"]] = row["EC_number"]
			prot_dict[row["EC_number"]].append(row["Parent_PDB_structure"])
		
		ec_num = []
		for i, row in df.iterrows():
			ec_num.append(ec_dict[row["Protein_file"].split("_")[0]])
	
		ecwise_count = Counter(ec_num)
		#print(ecwise_count)
		
		X_train_df = pd.DataFrame()
		Y_train_df = pd.DataFrame()
		X_test_df = pd.DataFrame()
		Y_test_df = pd.DataFrame()
		
		train_idx = []
		test_idx = []
		for ec, count in ecwise_count.items():
			prot_list = prot_dict[ec]
			sub_df = df[df["Protein_file"].str.contains("|".join(prot_list))]
			print(ec, len(sub_df.index), count)
			
			X_dataset = sub_df[list(df.columns)[:-1]]
			Y_dataset = sub_df[list(df.columns)[-1]]
			random_state = 42  #XXX: Very very important for reproducibility of the split
			
			if(len(sub_df.index)>=3):
				X_train, X_test, Y_train, Y_test = train_test_split(X_dataset, Y_dataset, test_size=0.1, random_state=random_state, shuffle=True)
				#print(len(X_train), len(X_test), len(Y_train), len(Y_test))
				train_idx.extend(list(X_train.index))
				test_idx.extend(list(X_test.index))
			elif(len(sub_df.index)==2):
				train_idx.append(list(sub_df.index)[0])
				test_idx.append(list(sub_df.index)[1])
			elif(len(sub_df.index)==1):
				train_idx.append(list(sub_df.index)[0])
		
		train_df = df.iloc[train_idx]
		test_df = df.iloc[test_idx]
		X_train_df = train_df[list(df.columns)[:-1]]
		X_test_df = test_df[list(df.columns)[:-1]]
		Y_train_df = train_df[list(df.columns)[-1]]
		Y_test_df = test_df[list(df.columns)[-1]]
		
		X_train_df.to_csv(outpath+act_type+"_X_train.csv", sep="\t", header=True, index=False)
		Y_train_df.to_csv(outpath+act_type+"_Y_train.csv", sep="\t", header=True, index=False)
		X_test_df.to_csv(outpath+act_type+"_X_test.csv", sep="\t", header=True, index=False)
		Y_test_df.to_csv(outpath+act_type+"_Y_test.csv", sep="\t", header=True, index=False)
				


























