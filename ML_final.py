import multiprocessing as mp
import sys

def main(args):
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

    tf = args[0]
    number = args[1]

    p = pd.read_csv(f"/u/scratch/s/seanchea/tfchip/cobindNew/qval/{tf}_train.csv", index_col = 0)
    def shuffle(i):
        row = p.iloc[i, :].values
        new = random.sample(list(row), len(row))
        df = pd.DataFrame()
        df[p.index[i]] = new
        return df

    a = []
    for i in range(len(p)):
        a.append(shuffle(i))
    n = pd.concat(a, 1)
    p = p.T
    p["y"] = [1] * len(p)
    n["y"] = [0] * len(p)
    alldata = pd.concat([p, n])
    y = alldata.y.tolist()
    X = alldata.drop("y", axis = 1)
    features = X.columns.tolist()
    X.columns = features
    X.shape[0] / 2

    p = pd.read_csv(f"/u/scratch/s/seanchea/tfchip/cobindNew/qval/{tf}_test.csv", index_col = 0)

    a = []
    for i in range(len(p)):
        a.append(shuffle(i))
    n = pd.concat(a, 1)
    p = p.T
    p["y"] = [1] * len(p)
    n["y"] = [0] * len(p)
    alldata = pd.concat([p, n])
    y_test = alldata.y.tolist()
    X_test = alldata.drop("y", axis = 1)
    X_test.columns = features

    normAnnoDir = f"/u/scratch/s/seanchea/tfchip/normalizedAnnoSALL2/qval/{tf}/"
    if tf=="SALL2":
        test = 0
        train = 1
    else:
        test = 1
        train = 0
    normAnnoBeds = os.listdir(normAnnoDir)
    normAnnoBeds.sort(key=lambda f: int(re.sub('\D', '', f)))
    annotations = pd.read_csv(normAnnoDir + normAnnoBeds[train], sep="\t", header = None)
    annotations[3] = [e[1] for e in annotations[3].str.split("peak")]
    annotations[3] = annotations[3].astype(int)
    annotations.sort_values(by=3, inplace=True)
    annotations.drop_duplicates(subset=3, inplace=True)
    distal_train = list(np.abs(annotations.iloc[:,7]) > 3000) * 2
    proximal_train = list(np.abs(annotations.iloc[:,7]) < 3000)* 2
    allpeaks_train = [True] * len(distal_train)
    if len(normAnnoBeds) > 1:
        annotations = pd.read_csv(normAnnoDir + normAnnoBeds[test], sep="\t", header = None)
        annotations[3] = [e[1] for e in annotations[3].str.split("peak")]
        annotations[3] = annotations[3].astype(int)
        annotations.sort_values(by=3, inplace=True)
        annotations.drop_duplicates(subset=3, inplace=True)
    distal_test = list(np.abs(annotations.iloc[:,7]) > 3000) * 2
    proximal_test = list(np.abs(annotations.iloc[:,7]) < 3000)* 2
    allpeaks_test = [True] * len(distal_test)

    contexts = [(distal_train, distal_test), (proximal_train, proximal_test), (allpeaks_train, allpeaks_test)]
    for index in range(len(contexts)):   
        X_context = X.loc[contexts[index][0], :]
        X_context_test = X_test.loc[contexts[index][1], :]
        y_context = list(compress(y, contexts[index][0]))
        y_context_test= list(compress(y_test, contexts[index][1]))

        param = {'force_row_wise': True, 'objective': 'binary', 'verbose': -1, 'metric': 'binary_logloss', 'num_leaves': 31}
        train_data = lgb.Dataset(X_context.iloc[0:len(X_context.index) - 1, ], label = y_context[0: len(y_context) - 1])
        validation_data = train_data.create_valid(X_context_test, label = y_context_test)
        start = time.time()
        trees = lgb.train(param, train_data, 10000, valid_sets = validation_data, early_stopping_rounds = 10, feature_name=features, verbose_eval=1)
        print(time.time() - start)

        print(f'Train set: {np.mean((trees.predict(X_context) > 0.5) == y_context)}')
        print(f'Test set: {np.mean((trees.predict(X_context_test) > 0.5) == y_context_test)}')
        print(f'Full train set: {np.mean((trees.predict(X) > 0.5) == y)}')
        print(f'Full test set: {np.mean((trees.predict(X_test) > 0.5) == y_test)}')

        start = time.time()
        explainer = shap.TreeExplainer(trees)
        shap_values_context = explainer.shap_values(X_context)
        print(f'Done with context {index} context-specific SHAP values in {time.time() - start} sec')
        if index != 2:
            start = time.time()
            shap_values_full = explainer.shap_values(X)
            print(f'Done with context {index} full SHAP values in {time.time() - start} sec')
        else:
            shap_values_full = shap_values_context

        cores = mp.cpu_count()
        chunks = np.array_split(X_context, cores)
        start = time.time()
        with concurrent.futures.ProcessPoolExecutor(cores) as executor:
            a = executor.map(explainer.shap_interaction_values, chunks)
        allInteractions_context = np.concatenate(list(a), axis = 0)
        print(f'Done with context {index} context-specific interactions in {time.time() - start} sec')

        chunks = np.array_split(X, cores)
        start = time.time()
        if index != 2:
            with concurrent.futures.ProcessPoolExecutor(cores) as executor:
                a = executor.map(explainer.shap_interaction_values, chunks)
            allInteractions_full = np.concatenate(list(a), axis = 0)
            print(f'Done with context {index} full interactions in {time.time() - start} sec')
        else:
            allInteractions_full = allInteractions_context

        saveDir=f'/u/scratch/s/seanchea/models/{tf}/{number}/python/context' + str(index)
        os.makedirs(saveDir, exist_ok = True)
        my_shelf = shelve.open(saveDir + "/saveData" ,'n')

        my_shelf["trees"] = locals()["trees"]
        my_shelf["X_context"] = locals()["X_context"]
        my_shelf["X_context_test"] = locals()["X_context_test"]
        my_shelf["allInteractions_context"] = locals()["allInteractions_context"]
        my_shelf["shap_values_context"] = locals()["shap_values_context"]
        my_shelf["X"] = locals()["X"]
        my_shelf["X_test"] = locals()["X_test"]
        my_shelf["allInteractions_full"] = locals()["allInteractions_full"]
        my_shelf["shap_values_full"] = locals()["shap_values_full"]
        my_shelf["features"] = locals()["features"]
        my_shelf["explainer"] = locals()["explainer"]

        my_shelf.close()

        print(f'Done saving all variables for context {index}')

if __name__ == '__main__':
    main(sys.argv[1:])
