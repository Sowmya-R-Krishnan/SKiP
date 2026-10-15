#Program to compute baseline performance metrics for the kinetics model

import sys
import csv
import pandas as pd
import numpy as np
import pickle
from sklearn import linear_model, model_selection, metrics, preprocessing, pipeline, neighbors, tree, ensemble, neural_network, kernel_ridge
from sklearn.tree import DecisionTreeRegressor
from sklearn.kernel_ridge import KernelRidge
from sklearn.neural_network import MLPRegressor
from sklearn.ensemble import GradientBoostingRegressor, RandomForestRegressor, AdaBoostRegressor
from sklearn.tree import DecisionTreeRegressor
from sklearn.neighbors import KNeighborsRegressor
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import mean_absolute_error, make_scorer, r2_score
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline
from scipy import stats
from scipy.stats import pearsonr

import catboost
from catboost import CatBoostRegressor

import xgboost
from xgboost import XGBRegressor

import lightgbm
from lightgbm import LGBMRegressor

inpath = sys.argv[1]
model_type = sys.argv[2]  #Km or kcat
nomit = sys.argv[3]
outfile = sys.argv[4]

#XXX: python 7_train_model_v1.py "./Sample_output/Random_split/" "kcat" 3 "./Sample_output/Random_split/kcat_LGBM_v1.pkl"

nomit = int(nomit)
X_train = pd.read_csv(inpath+model_type+"_X_train.csv", sep="\t", header=0)
X_test = pd.read_csv(inpath+model_type+"_X_test.csv", sep="\t", header=0)
Y_train = pd.read_csv(inpath+model_type+"_Y_train.csv", sep="\t", header=0)
Y_test = pd.read_csv(inpath+model_type+"_Y_test.csv", sep="\t", header=0)
print(len(X_train.index), len(Y_train.index), len(X_test.index), len(Y_test.index))

feats_final = list(X_train.columns)[nomit:]

def corr_function(Y1, Y2):
	corr, p_val = pearsonr(Y1, Y2)
	return corr


#TODO: Uncomment the model of interest and train it
def regression_model_without_norm(data, combination, y):
	X = data[combination]

	#model = linear_model.LinearRegression()
	#tuner = GridSearchCV(model, {}, scoring='r2', refit=True, cv=10)  #Since refit=True, the best model is automatically returned by the GridSearch function
	
	#model = linear_model.Ridge()
	#parameters = {'alpha':[1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0], 'solver':('auto', 'svd', 'cholesky', 'lsqr', 'sparse_cg', 'lbfgs')}
	#tuner = GridSearchCV(model, parameters, scoring=make_scorer(corr_function), refit=True, cv=10)
	
	#model = KNeighborsRegressor(n_jobs=20)
	#parameters = {'n_neighbors':[3, 4, 5, 6, 7, 8, 9, 10], 'weights':('uniform', 'distance')}
	#tuner = GridSearchCV(model, parameters, scoring=make_scorer(corr_function), refit=True, cv=10)
	
	#model = linear_model.HuberRegressor()
	#parameters = {'epsilon':[1.0, 1.15, 1.25, 1.35, 1.45], 'alpha':[0.0001, 0.001, 0.01, 0.1, 1.0], 'max_iter':[100, 200, 300]}
	#tuner = GridSearchCV(model, parameters, scoring=make_scorer(corr_function), refit=True, cv=10)
	
	#model = DecisionTreeRegressor(random_state=42)
	#parameters = {'criterion':('squared_error', 'absolute_error')}
	#tuner = GridSearchCV(model, parameters, scoring='r2', refit=True, cv=10)
	
	#model = RandomForestRegressor(random_state=42, oob_score=True, n_jobs=4)
	#parameters = {}
	#tuner = GridSearchCV(model, parameters, scoring='r2', refit=True, cv=10)
	
	#model = AdaBoostRegressor(random_state=42)
	#parameters = {'n_estimators':[50,100,150,200], 'learning_rate':[0.5,1.0,1.5,2.0], 'loss':('linear', 'square', 'exponential')}
	#tuner = GridSearchCV(model, parameters, scoring=make_scorer(corr_function), refit=True, cv=10)
	
	#model = GradientBoostingRegressor(random_state=42)
	#parameters = {'loss':('squared_error', 'absolute_error', 'huber'), 'n_estimators':[100, 200, 300, 400, 500], 'criterion':('friedman_mse', 'squared_error')}
	#tuner = GridSearchCV(model, parameters, scoring='r2', refit=True, cv=10)
	
	#model = CatBoostRegressor(loss_function='RMSE', random_state=42, verbose=1000)  #Write to stdout once in every 1000 epochs
	#parameters = {}
	#tuner = GridSearchCV(model, parameters, scoring='r2', refit=True, cv=10)
	
	#model = XGBRegressor(n_jobs=4, random_state=42)
	#parameters = {'n_estimators':[100,200,300], 'tree_method':('exact', 'approx', 'hist')}
	#tuner = GridSearchCV(model, parameters, scoring='r2', refit=True, cv=10)
	
	model = LGBMRegressor(random_state=42, n_jobs=4)
	parameters = {'boosting_type':['gbdt','dart','rf'], 'learning_rate':[0.1, 0.01, 0.001], 'n_estimators':[100,200,300,400,500]}
	tuner = GridSearchCV(model, parameters, scoring='r2', refit=True, cv=10)
	
	#model = KernelRidge()
	#parameters = {'kernel':['additive_chi2', 'chi2', 'linear', 'poly', 'polynomial', 'rbf', 'laplacian', 'sigmoid', 'cosine']}
	#tuner = GridSearchCV(model, parameters, scoring='r2', refit=True, cv=10)

	#model = MLPRegressor(random_state=42)
	#parameters = {'hidden_layer_sizes':[100,200,300,400,500], 'activation':('tanh', 'relu'), 'batch_size':[8,32,64,128,256], 'learning_rate':('constant', 'adaptive'), 'learning_rate_init':[0.001, 0.01, 0.005], 'max_iter':[200,400,600,800,1000]}
	#tuner = GridSearchCV(model, parameters, scoring=make_scorer(corr_function), refit=True, cv=10)
	
	tuner.fit(X, y)
	predicted = tuner.predict(X)
	r2 = r2_score(y, predicted)
	pearson_corr, p_value = pearsonr(y, predicted)
	mae = mean_absolute_error(y, predicted)
	return(r2, pearson_corr, mae, tuner)
	
train_r2, train_corr, train_mae, model = regression_model_without_norm(X_train, feats_final, Y_train[model_type+"_value"])

test_preds = model.predict(X_test[feats_final])
test_r2 = r2_score(Y_test[model_type+"_value"], test_preds)
test_corr, _ = pearsonr(Y_test[model_type+"_value"], test_preds)
test_mae = mean_absolute_error(Y_test[model_type+"_value"], test_preds)
	
with open(outfile, "wb") as out:
	pickle.dump(model, file=out)
out.close()

print("Train R2: "+str(np.round(train_r2, 5))+", Train r: "+str(np.round(train_corr, 5))+", Train MAE: "+str(np.round(train_mae, 5))+", Test R2: "+str(np.round(test_r2, 5))+", Test r: "+str(np.round(test_corr, 5))+", Test MAE: "+str(np.round(test_mae, 5)))


























