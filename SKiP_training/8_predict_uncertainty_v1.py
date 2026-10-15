#Program to extract aleatoric and epistemic uncertainty for a given dataset

import os
import sys
import csv
import pickle
import pandas as pd
import numpy as np
from mapie import regression
from mapie.regression import MapieRegressor
from mapie.metrics import regression_coverage_score

infile = sys.argv[1]
nomit = sys.argv[2]
model_type = sys.argv[3]
cv_models = sys.argv[4]
mapie_results = sys.argv[5]
outfile = sys.argv[6]
gt = sys.argv[7]

#XXX: python 8_predict_uncertainty_v1.py "./Sample_output/Random_split/kcat_X_test.csv" 3 "kcat" "./Final_SKiP_models/kcat_LGBM_10F_CV/" "./Final_SKiP_models/kcat_LGBM_MAPIE_results_v1.pkl" "./Sample_output/kcat_X_test_uncertainty_sample.csv" "./Sample_output/Random_split/kcat_Y_test.csv"

nomit = int(nomit)
data = pd.read_csv(infile, sep="\t", header=0)
data_y = pd.read_csv(gt, sep="\t", header=0)
print(len(data.index))

feats_final = list(data.columns)[nomit:]
def corr_function(Y1, Y2):
	corr, p_val = pearsonr(Y1, Y2)
	return corr
	
with open(mapie_results, "rb") as f:
	mapie = pickle.load(f)
	
#Compute epistemic uncertainty with the 10-fold CV models
cv_preds = [[] for _ in range(len(data.index))]
for model_file in os.listdir(cv_models):
	with open(cv_models+model_file, "rb") as f:
		model = pickle.load(f)
		
	preds = model.predict(data[feats_final])
	
	for i, row in data.iterrows():
		cv_preds[i].append(preds[i])
		
epistemic_error = [np.var(l) for l in cv_preds]  #XXX: IMPORTANT NOTE: Only Variance should be used, not standard deviation

#Compute aleatoric uncertainty with MAPIE results
pred_mean, pred_intervals = mapie.predict(data[feats_final], alpha=0.1)

with open(outfile, "w") as out:
	print("Protein_file\tLigand_file\tES_complex\tActual_p"+model_type+"\tEpistemic_error\tMAPIE_mean_pred\tMAPIE_CI_lower\tMAPIE_CI_upper", file=out)
	for i, row in data.iterrows():
		print(str(row["Protein_file"])+"\t"+str(row["name"])+"\t"+str(row["ES_complex"])+"\t"+str(data_y.iloc[i][model_type+"_value"])+"\t"+str(epistemic_error[i])+"\t"+str(pred_mean[i])+"\t"+str(pred_intervals[i][0])+"\t"+str(pred_intervals[i][1]), file=out)
out.close()

coverage_score = regression_coverage_score(data_y[model_type+"_value"], pred_intervals[:, 0, 0], pred_intervals[:, 1, 0])
print("Regression coverage score: "+str(coverage_score))


























