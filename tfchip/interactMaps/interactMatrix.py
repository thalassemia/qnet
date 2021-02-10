import multiprocessing as mp
mp.set_start_method('spawn', force = True)
import shelve
import pandas as pd
import numpy as np
import shap
import concurrent.futures
from itertools import compress
import seaborn as sns
import matplotlib
matplotlib.use('SVG')
import matplotlib.pylab as plt
import scipy as sc
import os

def load_data(dirname):
    filename='/u/scratch/s/seanchea/models/' + dirname + '/' + str(0) + '/python/self/saveData'
    my_shelf = shelve.open(filename)
    for key in my_shelf:
        globals()[key]=my_shelf[key]
    my_shelf.close()
    pos = int(allInteractions.shape[0]/2)
    avg = [np.mean(allInteractions[:pos, :, :], axis=0)]
    maximum = allInteractions[:pos, :, :]/100
    for i in range(99):
        filename='/u/scratch/s/seanchea/models/' + dirname + '/' + str(i + 1) + '/python/self/saveData'
        my_shelf = shelve.open(filename)
        for key in my_shelf:
            globals()[key]=my_shelf[key]
        my_shelf.close()
        avg1 = np.mean(allInteractions[:pos, :, :], axis=0)
        avg.append(avg1)
        maximum = maximum + allInteractions[:pos, :, :]/100
        print(i + 1)
    avg = np.mean(np.stack(avg, axis = 0), axis = 0)
    avgMax = np.amax(maximum, axis = 0)
    return [avg, avgMax, features]

os.makedirs("/u/home/s/seanchea/models/plots/all_train.csv/self/", exist_ok = True)  
data = load_data("all_train.csv")
np.savetxt("/u/home/s/seanchea/plots/all_train.csv/self/max.txt", data[1], delimiter=",")
np.savetxt("/u/home/s/seanchea/plots/all_train.csv/self/avg.txt", data[0], delimiter=",")