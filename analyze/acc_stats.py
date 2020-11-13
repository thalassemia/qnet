import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

train_chea = np.zeros(100)
test_chea = np.zeros(100)

with open("/u/home/s/seanchea/ml_logs/chea.log") as f:
    line = f.readline()
    i = 0
    while line:
        if "Train set: " in line:
            train_chea[i] = line.strip().split("Train set: ")[1]
        elif "Test set: " in line:
            test_chea[i] = line.strip().split("Test set: ")[1]
            i = i + 1
        line = f.readline()

train_lisa = np.zeros(100)
test_lisa = np.zeros(100)

with open("/u/home/s/seanchea/ml_logs/lisa.log") as f:
    line = f.readline()
    i = 0
    while line:
        if "Train set: " in line:
            train_lisa[i] = line.strip().split("Train set: ")[1]
        elif "Test set: " in line:
            test_lisa[i] = line.strip().split("Test set: ")[1]
            i = i + 1
        line = f.readline()

train_up = np.zeros(100)
test_up = np.zeros(100)

with open("/u/home/s/seanchea/ml_logs/uropa_up.log") as f:
    line = f.readline()
    i = 0
    while line:
        if "Train set: " in line:
            train_up[i] = line.strip().split("Train set: ")[1]
        elif "Test set: " in line:
            test_up[i] = line.strip().split("Test set: ")[1]
            i = i + 1
        line = f.readline()

accs = [np.mean(train_chea), np.mean(test_chea), np.mean(train_lisa), np.mean(test_lisa), np.mean(train_up), np.mean(test_up)]
error = [np.std(train_chea), np.std(test_chea), np.std(train_lisa), np.std(test_lisa), np.std(train_up), np.std(test_up)]

np.savetxt('/u/home/s/seanchea/ml_logs/accs.txt', accs, delimiter = ',')
np.savetxt('/u/home/s/seanchea/ml_logs/errs.txt', error, delimiter = ',')