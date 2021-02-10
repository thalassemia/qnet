import multiprocessing as mp
import sys

def main(args):
    # this flag is needed because the default "fork" does not work here
    mp.set_start_method("spawn")
    import concurrent.futures
    import pandas as pd
    import numpy as np
    import time
    import random
    from itertools import compress
    import os
    import shap
    import lightgbm as lgb
    import shelve
    import re

    def model(number, exp, control):
        print(exp)
        p = pd.read_csv("/u/scratch/s/seanchea/cobinding/cobind/" + exp, sep = "\t", index_col = 0)
        # independently shuffle each column, theoretically eliminating any TF binding dependencies
        n = p.apply(lambda x: x.sample(frac = 1).values) # .values needed or else will return Series that realigns via index

        p["y"] = [1] * len(p)
        n["y"] = [0] * len(p)
        alldata = pd.concat([p, n])
        y = alldata.y.tolist()
        X = alldata.drop("y", axis = 1)
        features = X.columns.tolist()
        X.columns = features
        X.shape[0] / 2

        p = pd.read_csv("/u/scratch/s/seanchea/cobinding/cobind/" + control, sep = "\t", index_col = 0)
        n = p.apply(lambda x: x.sample(frac = 1).values)

        p["y"] = [1] * len(p)
        n["y"] = [0] * len(p)
        alldata = pd.concat([p, n])
        y_test = alldata.y.tolist()
        X_test = alldata.drop("y", axis = 1)
        X_test.columns = features

        param = {'force_row_wise': True, 'objective': 'binary', 'verbose': -1, 'metric': 'binary_logloss', 'num_leaves': 31}
        # lightGBM does not work when the number of positive and negative sample is exactly equal so arbitarily remove final negative sample
        train_data = lgb.Dataset(X.iloc[0:len(X.index) - 1, ], label = y[0: len(y) - 1])
        validation_data = train_data.create_valid(X_test, label = y_test)
        start = time.time()
        trees = lgb.train(param, train_data, 10000, valid_sets = validation_data, early_stopping_rounds = 50, feature_name=features, verbose_eval=1)
        print(f'Done training in {time.time() - start} sec')

        #TO-DO: implement more sophisticated evaluation measure like AUROC
        print(f'Train set: {np.mean((trees.predict(X) > 0.5) == y)}')
        print(f'Test set: {np.mean((trees.predict(X_test) > 0.5) == y_test)}')

        start = time.time()
        explainer = shap.TreeExplainer(trees)
        shap_values = explainer.shap_values(X)
        print(f'Done with SHAP values in {time.time() - start} sec')

        cores = mp.cpu_count()
        chunks = np.array_split(X, cores)
        start = time.time()
        # having SHAP calculate interaction values on chunks to enable multiprocessing
        with concurrent.futures.ProcessPoolExecutor(cores) as executor:
            a = executor.map(explainer.shap_interaction_values, chunks)
        allInteractions = np.concatenate(list(a), axis = 0)
        print(f'Done with interactions in {time.time() - start} sec')

        # export all relevant variables for future reference
        saveDir=f'/u/scratch/s/seanchea/models/{exp}/{number}/python/test/'
        os.makedirs(saveDir, exist_ok = True)
        my_shelf = shelve.open(saveDir + "/saveData" ,'n', protocol = 4)

        my_shelf["trees"] = locals()["trees"]
        my_shelf["X"] = locals()["X"]
        my_shelf["X_test"] = locals()["X_test"]
        my_shelf["allInteractions"] = locals()["allInteractions"]
        my_shelf["shap_values"] = locals()["shap_values"]
        my_shelf["features"] = locals()["features"]
        my_shelf["explainer"] = locals()["explainer"]

        my_shelf.close()

        print('Done saving all variables (with test set)')

        # repeat procedure for model validated on itself
        start = time.time()
        trees = lgb.train(param, train_data, 100, valid_sets = train_data, feature_name=features, verbose_eval=1)
        print(f'Done training in {time.time() - start} sec')

        #TO-DO: implement more sophisticated evaluation measure like AUROC
        print(f'Train set: {np.mean((trees.predict(X) > 0.5) == y)}')
        print(f'Test set: {np.mean((trees.predict(X_test) > 0.5) == y_test)}')

        start = time.time()
        explainer = shap.TreeExplainer(trees)
        shap_values = explainer.shap_values(X)
        print(f'Done with SHAP values in {time.time() - start} sec')

        cores = mp.cpu_count()
        chunks = np.array_split(X, cores)
        start = time.time()
        # having SHAP calculate interaction values on chunks to enable multiprocessing
        with concurrent.futures.ProcessPoolExecutor(cores) as executor:
            a = executor.map(explainer.shap_interaction_values, chunks)
        allInteractions = np.concatenate(list(a), axis = 0)
        print(f'Done with interactions in {time.time() - start} sec')

        # export all relevant variables for future reference
        saveDir=f'/u/scratch/s/seanchea/models/{exp}/{number}/python/self/'
        os.makedirs(saveDir, exist_ok = True)
        my_shelf = shelve.open(saveDir + "/saveData" ,'n', protocol = 4)

        my_shelf["trees"] = locals()["trees"]
        my_shelf["X"] = locals()["X"]
        my_shelf["X_test"] = locals()["X_test"]
        my_shelf["allInteractions"] = locals()["allInteractions"]
        my_shelf["shap_values"] = locals()["shap_values"]
        my_shelf["features"] = locals()["features"]
        my_shelf["explainer"] = locals()["explainer"]

        my_shelf.close()

        print('Done saving all variables')

    for i in range(100): 
        model(i, args[-2], args[-1])

# this convoluted setup is necessary to enable the "spawn" flag to be set properly
if __name__ == '__main__':
    main(sys.argv[:])