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

anno_dir = os.path.expandvars("$SCRATCH/tfchip/cheaannotated/")
normDir = os.path.expandvars("$SCRATCH/tfchip/Cnormalized/")
intDir = os.path.expandvars("$SCRATCH/tfchip/Cintersect/")
cobindDir = os.path.expandvars("$SCRATCH/tfchip/Ccobind/")
deFile = os.path.expandvars("$SCRATCH/data/bsf.csv")
de = pd.read_csv(deFile, sep = "\t", header = None)
de = de.loc[de[2] < 0.05]
de.rename(columns = {0:"gene_name"}, inplace = True)
de.drop_duplicates(subset = "gene_name", inplace = True)

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

with concurrent.futures.ProcessPoolExecutor(36) as executor:
    beds = glob.glob(anno_dir + "*/*.csv")
    list(tqdm(executor.map(norm, beds), total=len(beds)))

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

with concurrent.futures.ProcessPoolExecutor(36) as executor:
    beds = glob.glob(normDir + "*/*.csv")
    list(tqdm(executor.map(intersect, beds), total = len(beds)))

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

cobind = pd.DataFrame()
deGenes = []

for i in range(len(sortedBeds)):
    factor = sortedBeds[i][0].split("/")[-2]
    temp = pd.read_csv(sortedBeds[i][0], sep = "\t")
    temp.drop_duplicates(subset = "gene_name", inplace = True, keep = "first", ignore_index = True)
    deGenes = temp["gene_name"].to_list()
    cobind[factor] = temp["peak_score"]
cobind.index = deGenes
os.makedirs(cobindDir, exist_ok = True)
cobind.to_csv(os.path.join(cobindDir, "chea_train.csv"), sep = "\t")

cobind = pd.DataFrame()
deGenes = []
for i in range(len(sortedBeds)):
    factor = sortedBeds[i][0].split("/")[-2]
    temp = pd.read_csv(sortedBeds[i][0], sep = "\t")
    if len(sortedBeds[i]) > 1:
        temp = pd.read_csv(sortedBeds[i][1], sep = "\t")
    temp.drop_duplicates(subset = "gene_name", inplace = True, keep = "first", ignore_index = True)
    deGenes = temp["gene_name"].to_list()
    cobind[factor] = temp["peak_score"]
cobind.index = deGenes
cobind.to_csv(os.path.join(cobindDir, "chea_test.csv"), sep = "\t")

anno_dir = os.path.expandvars("$SCRATCH/tfchip/lisaannotated/")
normDir = os.path.expandvars("$SCRATCH/tfchip/Lnormalized/")
intDir = os.path.expandvars("$SCRATCH/tfchip/Lintersect/")
cobindDir = os.path.expandvars("$SCRATCH/tfchip/Lcobind/")

with concurrent.futures.ProcessPoolExecutor(36) as executor:
    beds = glob.glob(anno_dir + "*/*.csv")
    list(tqdm(executor.map(norm, beds), total=len(beds)))

with concurrent.futures.ProcessPoolExecutor(36) as executor:
    beds = glob.glob(normDir + "*/*.csv")
    list(tqdm(executor.map(intersect, beds), total = len(beds)))

factors = os.listdir(intDir)
sortedBeds = []
for factor in factors:
    beds = glob.glob(os.path.join(intDir, factor, "*.csv"))
    beds.sort(key = num)
    sortedBeds.append(beds)

cobind = pd.DataFrame()
deGenes = []

for i in range(len(sortedBeds)):
    factor = sortedBeds[i][0].split("/")[-2]
    temp = pd.read_csv(sortedBeds[i][0], sep = "\t")
    temp.drop_duplicates(subset = "gene_name", inplace = True, keep = "first", ignore_index = True)
    deGenes = temp["gene_name"].to_list()
    cobind[factor] = temp["peak_score"]
cobind.index = deGenes
os.makedirs(cobindDir, exist_ok = True)
cobind.to_csv(os.path.join(cobindDir, "lisa_train.csv"), sep = "\t")

cobind = pd.DataFrame()
deGenes = []
for i in range(len(sortedBeds)):
    factor = sortedBeds[i][0].split("/")[-2]
    temp = pd.read_csv(sortedBeds[i][0], sep = "\t")
    if len(sortedBeds[i]) > 1:
        temp = pd.read_csv(sortedBeds[i][1], sep = "\t")
    temp.drop_duplicates(subset = "gene_name", inplace = True, keep = "first", ignore_index = True)
    deGenes = temp["gene_name"].to_list()
    cobind[factor] = temp["peak_score"]
cobind.index = deGenes
cobind.to_csv(os.path.join(cobindDir, "lisa_test.csv"), sep = "\t")