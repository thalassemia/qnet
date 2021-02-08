import multiprocessing as mp
mp.set_start_method('spawn', force = True)
import shelve
import pandas as pd
import numpy as np
import shap
import concurrent.futures
from itertools import compress
import seaborn as sns
import matplotlib.pylab as plt
import scipy as sc
import os

def load_data(dirname):
    filename='/u/scratch/s/seanchea/atac/' + dirname + '/' + str(0) + '/python/saveData'
    my_shelf = shelve.open(filename)
    for key in my_shelf:
        globals()[key]=my_shelf[key]
    my_shelf.close()
    pos = int(allInteractions.shape[0]/2)
    avg = [np.mean(allInteractions[:pos, :, :], axis=0)]
    maximum = allInteractions[:pos, :, :]/100
    for i in range(100):
        filename='/u/scratch/s/seanchea/atac/' + dirname + '/' + str(i + 1) + '/python/saveData'
        my_shelf = shelve.open(filename)
        for key in my_shelf:
            globals()[key]=my_shelf[key]
        my_shelf.close()
        avg1 = np.mean(allInteractions[:pos, :, :], axis=0)
        avg.append(avg1)
        maximum = maximum + allInteractions[:pos, :, :]/20
        print(i + 1)
    avg = np.mean(np.stack(avg, axis = 0), axis = 0)
    avgMax = np.amax(maximum, axis = 0)
    return [avg, avgMax, features]

def heatmap(features, mat, threshold, outname):
    sigFeatures = []
    include = [False] * len(features)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            if mat[i,j] > threshold:
                include[i] = True
                sigFeatures.append(features[i])
                break

    pruned = mat[include, :]
    pruned = pruned[:, include]
    sns.set(font_scale=1)
    ax = sns.clustermap(pruned, method="ward", center=0, xticklabels=sigFeatures, yticklabels=sigFeatures, cbar_kws={'orientation': 'horizontal'}, dendrogram_ratio=0.05, cbar_pos = (0.20, 0, 0.5,0.02))
    ax.ax_heatmap.set_xticklabels(ax.ax_heatmap.get_xmajorticklabels(), fontsize = 14, rotation = 45)
    ax.ax_heatmap.set_yticklabels(ax.ax_heatmap.get_ymajorticklabels(), fontsize = 14)
    ax.ax_heatmap.tick_params(right=False, bottom=False)
    ax.savefig(f"/u/scratch/s/seanchea/atac/plots/{outname}.svg", dpi = 1200)

    iVals = []
    name = []
    for i in range(mat.shape[0]):
        for j in range(i, mat.shape[1]):
            iVals.append(mat[i,j])
            if i != j:
                name.append(features[i] + "-" + features[j])
            else:
                name.append(features[i])
    if outname.split("_")[1] == "max":
        unrolled = pd.DataFrame({"Max Interaction Value": iVals})
        unrolled.index = name
        top = unrolled.sort_values(by="Max Interaction Value")[-10:]
    else:
        unrolled = pd.DataFrame({"Avg. Interaction Value": iVals})
        unrolled.index = name
        top = unrolled.sort_values(by="Avg. Interaction Value")[-10:]
    ax = top.plot.barh()
    fig = ax.get_figure()
    plt.tight_layout()
    fig.savefig(f"/u/scratch/s/seanchea/atac/plots/{outname}_bar.svg", dpi = 1200)

os.makedirs("/u/scratch/s/seanchea/atac/plots/", exist_ok = True)  
data = load_data("up_models")
np.savetxt("/u/home/s/seanchea/plots/umax.txt", data[1], delimiter=",")
np.savetxt("/u/home/s/seanchea/plots/uavg.txt", data[0], delimiter=",")
heatmap(data[2], data[1], 0.1, "up_max")
heatmap(data[2], data[0], 0.01, "up_avg")