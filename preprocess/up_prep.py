import multiprocessing as mp
import pandas as pd
import os
import subprocess
import json
import shlex
import concurrent.futures
from tqdm import tqdm
import glob
import numpy as np

anno_dir = os.path.expandvars("$SCRATCH/tfchip/annotated/")
normDir = os.path.expandvars("$SCRATCH/tfchip/Unormalized/")
intDir = os.path.expandvars("$SCRATCH/tfchip/Uintersect/")
cobindDir = os.path.expandvars("$SCRATCH/tfchip/Ucobind/")
deFile = os.path.expandvars("$SCRATCH/data/bsf.csv")
de = pd.read_csv(deFile, sep = "\t", header = None)
de = de.loc[de[2] < 0.05]
de.rename(columns = {0:"gene_name"}, inplace = True)
de.drop_duplicates(subset = "gene_name", inplace = True)
up = de.loc[de[1] > 0]
down = de.loc[de[1] < 0]

def norm(a):
    temp = pd.read_csv(a)
    try:
        temp.sort_values(by="V9", ascending=False, inplace=True)
    except:
        temp.sort_values(by="V5", ascending=False, inplace=True)
    temp.index = range(0,len(temp.index))
    temp["peak_score"] = (-1*(temp.index)+len(temp.index))/len(temp.index)
    temp.dropna(subset = ["SYMBOL"], inplace = True)
    bed = a.split("/")[-1]
    factor = a.split("/")[-2]
    os.makedirs(os.path.join(normDir, factor), exist_ok=True)
    temp.to_csv(os.path.join(normDir, factor, bed), sep="\t", index=False)

""" with concurrent.futures.ProcessPoolExecutor(36) as executor:
    beds = glob.glob(anno_dir + "*/*.csv")
    list(tqdm(executor.map(norm, beds), total=len(beds))) """

def intersect(a):
    bed = pd.read_csv(a, sep = "\t")
    bed.rename(columns = {"SYMBOL": "gene_name"}, inplace = True)
    merged = de.merge(bed, on = "gene_name", how = "left")
    merged = merged.loc[:, ["gene_name", "peak_score"]]
    merged.fillna(0, inplace = True)
    bedname = a.split("/")[-1]
    factor = a.split("/")[-2]
    os.makedirs(os.path.join(intDir, factor), exist_ok = True)
    merged.to_csv(os.path.join(intDir, factor, bedname), sep = "\t", index = False)

""" with concurrent.futures.ProcessPoolExecutor(36) as executor:
    beds = glob.glob(normDir + "*/*.csv")
    list(tqdm(executor.map(intersect, beds), total = len(beds))) """

def num(x):
    try:
        return int(x.split("/")[-1].split(".bed")[0].split("_")[0])
    except:
        return x

factors = os.listdir(intDir)
sortedBeds = []
for factor in factors:
    beds = glob.glob(os.path.join(intDir, factor, "*.csv"))
    beds.sort(key = num)
    sortedBeds.append(beds)

cobindUp = pd.DataFrame()
cobindDown = pd.DataFrame()
cobindAll = pd.DataFrame()
deGenes = []

for i in range(len(sortedBeds)):
    factor = sortedBeds[i][0].split("/")[-2]
    temp = pd.read_csv(sortedBeds[i][0], sep = "\t")
    temp.drop_duplicates(subset = "gene_name", inplace = True, keep = "first", ignore_index = True)
    if factor in up["gene_name"].to_list():
        cobindUp[factor] = temp["peak_score"]
    else:
        cobindDown[factor] = temp["peak_score"]
    deGenes = temp["gene_name"].to_list()
    cobindAll[factor] = temp["peak_score"]
cobindUp.index = deGenes
cobindDown.index = deGenes
cobindAll.index = deGenes
os.makedirs(cobindDir, exist_ok = True)
cobindUp.to_csv(os.path.join(cobindDir, "up_train.csv"), sep = "\t")
cobindDown.to_csv(os.path.join(cobindDir, "down_train.csv"), sep = "\t")
cobindAll.to_csv(os.path.join(cobindDir, "all_train.csv"), sep = "\t")


cobindUp = pd.DataFrame()
cobindDown = pd.DataFrame()
cobindAll = pd.DataFrame()
deGenes = []
for i in range(len(sortedBeds)):
    factor = sortedBeds[i][0].split("/")[-2]
    temp = pd.read_csv(sortedBeds[i][0], sep = "\t")
    if len(sortedBeds[i]) > 1:
        temp = pd.read_csv(sortedBeds[i][1], sep = "\t")
    temp.drop_duplicates(subset = "gene_name", inplace = True, keep = "first", ignore_index = True)
    if factor in up["gene_name"].to_list():
        cobindUp[factor] = temp["peak_score"]
    else:
        cobindDown[factor] = temp["peak_score"]
    deGenes = temp["gene_name"].to_list()
    cobindAll[factor] = temp["peak_score"]
cobindUp.index = deGenes
cobindDown.index = deGenes
cobindAll.index = deGenes
cobindUp.to_csv(os.path.join(cobindDir, "up_test.csv"), sep = "\t")
cobindDown.to_csv(os.path.join(cobindDir, "down_test.csv"), sep = "\t")
cobindAll.to_csv(os.path.join(cobindDir, "all_test.csv"), sep = "\t")