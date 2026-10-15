#Program to compute all performance metrics for a model

import os
import sys
import csv
import pandas as pd
import numpy as np
from scipy.stats import pearsonr, spearmanr
from sklearn.metrics import r2_score, mean_absolute_error

gt = sys.argv[1]
preds = sys.argv[2]
gt_col = sys.argv[3]
pred_col = sys.argv[4]

#XXX: python 15_compute_metrics_v1.py "./data/kcat_independent_test_set_v1.csv" "./Sample_output/SKiP_kcat_preds_v1.csv" "pkcat_value" "Pred_pkcat"

gt_df = pd.read_csv(gt, sep="\t", header=0)
pred_df = pd.read_csv(preds, sep="\t", header=0)

l1 = list(gt_df[gt_col])
l2 = list(pred_df[pred_col])

r2 = r2_score(l1, l2)
r = pearsonr(l1, l2)[0]
mae = mean_absolute_error(l1, l2)
rho = spearmanr(l1, l2)[0]

ae = [abs(v1-v2) for v1,v2 in zip(l1, l2)]
median_ae = np.median(ae)

print(f'R2 score = {r2}')
print(f'Pearson\'s r = {r}')
print(f'MAE = {mae}')
print(f'Median AE = {median_ae}')
print(f'Spearman\'s r = {rho}')























