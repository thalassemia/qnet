import multiprocessing as mp
mp.set_start_method("spawn", force = True)
import pandas as pd
import numpy as np
import concurrent.futures
from itertools import compress
import seaborn as sns
import matplotlib
matplotlib.use('SVG')
import matplotlib.pylab as plt
import scipy as sc
import os
import glob
import shelve
from tqdm import tqdm    

def heatmap(features, mat, outname, threshold):
    sigFeatures = []
    include = [False] * len(features)
    mat = np.loadtxt(mat, delimiter = ',')
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            if mat[i,j] > threshold:
                include[i] = True
                sigFeatures.append(features[i])
                break

    pruned = mat[include, :]
    pruned = pruned[:, include]
    sns.set(font_scale=1)
    ax = sns.clustermap(pruned, method="ward", center=0, xticklabels=sigFeatures, yticklabels=sigFeatures, cbar_kws={'orientation': 'horizontal'}, dendrogram_ratio=0.05, cbar_pos = (0.20, -0.05, 0.5,0.02))
    ax.ax_heatmap.tick_params(right=False, bottom=False)
    ax.ax_heatmap.set_xticklabels(ax.ax_heatmap.get_xmajorticklabels(), fontsize = 10, rotation = 90)
    ax.ax_heatmap.set_yticklabels(ax.ax_heatmap.get_ymajorticklabels(), fontsize = 10, rotation = 0)
    ax.savefig(f"/u/home/s/seanchea/models/plots/{outname}.svg", dpi = 1200)

    iVals = []
    name = []
    for i in range(mat.shape[0]):
        for j in range(i, mat.shape[1]):
            iVals.append(mat[i,j])
            if i != j:
                name.append(features[i] + "-" + features[j])
            else:
                name.append(features[i])
    if outname.split("/")[-1] == "max":
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
    fig.savefig(f"/u/home/s/seanchea/models/plots/{outname}_bar.svg", dpi = 1200)

def loadFeatures(mod):
    name, method = mod
    filename='/u/scratch/s/seanchea/models/' + name + '/' + str(0) + '/python/' + method + '/saveData'
    my_shelf = shelve.open(filename)
    features = my_shelf["features"]
    my_shelf.close()
    return features

def unzip(args):
    heatmap(*args)

if __name__ == '__main__':
    matrices = glob.glob("/u/home/s/seanchea/plots/**/*.txt", recursive = True)
    outnames = [i.split('/')[-3] + '/' + i.split('/')[-2] + '/' + i.split('/')[-1].split('.')[-2] for i in matrices]
    threshold = [0.1 if i.split('/')[-1].split('.')[-2] == "max" else 0.01 for i in matrices]

    feats = []
    model = [(i.split('/')[-3], i.split('/')[-2]) for i in matrices]
    with concurrent.futures.ProcessPoolExecutor(36) as executor:
        feats = list(tqdm(executor.map(loadFeatures, model), total = len(model)))
        list(tqdm(executor.map(unzip, zip(feats, matrices, outnames, threshold)), total = len(model)))