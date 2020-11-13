# Train/Interpret ML Models

**Dependencies:** Python 3.8, pandas 1.1.3, numpy 1.19.2, shap 0.36.0, lightgbm 3.0.0

Trains 100 LightGBM models on train and test sets with the following settings:

    'objective': 'binary'
    'metric': 'binary_logloss'
    early_stopping_rounds = 50

Calculates train and test set accuracies for each model (saved in Hoffman job output file).

Saves the LightGBM model, the training set, the test set, the SHAP interaction values, the SHAP values, the feature list (TF list), and the SHAP explainer object in a python shelf for future access.