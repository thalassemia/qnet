import pandas as pd
import numpy as np
import random
import os

p = pd.read_csv("/u/scratch/s/seanchea/tfchip/Ucobind/up_train.csv", sep = "\t", index_col = 0)
def shuffle(i):
    row = p.iloc[:, i].values
    new = random.sample(list(row), len(row))
    df = pd.DataFrame()
    df[p.columns[i]] = new
    return df

a = []
for i in range(p.shape[1]):
    a.append(shuffle(i))
n = pd.concat(a, 1)

n = n.T
n.to_csv("/u/scratch/s/seanchea/tfchip/Ucobind/up_train_negative.csv", sep = "\t")