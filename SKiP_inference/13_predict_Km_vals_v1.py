#Program to load the models and predict Km values

import sys
import csv
import pickle
import pandas as pd
import numpy as np
from sklearn.metrics import mean_absolute_error, make_scorer, r2_score
from scipy.stats import pearsonr

infile = sys.argv[1]
model = sys.argv[2]
nomit = sys.argv[3]
outfile = sys.argv[4]

#XXX: python 13_predict_Km_vals_v1.py "./Sample_output/Features/Km_independent_test_set_final_feats_v1.csv" "../SKiP_training/Final_SKiP_models/Km_RF_v1.pkl" 3 "./Sample_output/SKiP_Km_preds_v1.csv"

test_df = pd.read_csv(infile, sep="\t", header=0)

def corr_function(Y1, Y2):
	corr, p_val = pearsonr(Y1, Y2)
	return corr

nomit = int(nomit)
feat_cols = []
feat_cols.extend(list(test_df.columns)[nomit:])
test_X = test_df[feat_cols]

with open(model, "rb") as f:
	weights = pickle.load(f)

test_preds = weights.predict(test_X)

with open(outfile, "w") as out:
	print("Protein_file\tname\tES_complex\tPred_pKm", file=out)
	nomit = int(nomit)
	test_X = test_df[list(test_df.columns)[nomit:]]
	
	for i, row in test_df.iterrows():
			print(str(row["Protein_file"])+"\t"+str(row["name"])+"\t"+str(row["ES_complex"])+"\t"+str(test_preds[i]), file=out)

out.close()
	





























