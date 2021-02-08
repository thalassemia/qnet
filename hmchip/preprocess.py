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
import io

anno_dir = os.path.expandvars("$SCRATCH/hmchip/annotated/")
normDir = os.path.expandvars("$SCRATCH/hmchip/normalized/")
intDir = os.path.expandvars("$SCRATCH/hmchip/intersect/")
cobindDir = os.path.expandvars("$SCRATCH/hmchip/cobind/")
bedDir = os.path.expandvars("$SCRATCH/hmchip/cleanbeds/")
cleanPath = os.path.expandvars("$SCRATCH/atac/corr/atacAnno_clean.tsv")
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
    # Merge with DE gene list
    """ bed = pd.read_csv(a, sep = "\t")
    bed.rename(columns = {"SYMBOL": "gene_name"}, inplace = True)
    merged = de.merge(bed, on = "gene_name", how = "left")
    merged = merged.loc[:, ["gene_name", "peak_score"]]
    merged.fillna(0, inplace = True)
    bedname = a.split("/")[-1]
    factor = a.split("/")[-2]
    os.makedirs(os.path.join(intDir, factor), exist_ok = True)
    merged.to_csv(os.path.join(intDir, factor, bedname), sep = "\t", index = False) """

    # Merge with ATAC gene list
    bedname = a.split("/")[-1]
    factor = a.split("/")[-2]
    bed = pd.read_csv(a, sep = "\t")
    bed = bed.loc[:, ["seqnames", "start", "end", "peak_score"]]
    os.makedirs(os.path.join(bedDir, factor), exist_ok = True)
    bed.to_csv(os.path.join(bedDir, factor, bedname), sep = "\t", index = False, header = False)
    command = shlex.split(f"bedtools intersect -a \"{cleanPath}\" -b \"{os.path.join(bedDir, factor, bedname)}\" -loj")
    result = subprocess.run(command, stdout=subprocess.PIPE)
    output = result.stdout
    df = pd.read_csv(io.BytesIO(output), sep="\t", header = None)
    scoreCol = df.shape[1] - 1
    df.loc[df[scoreCol] == ".", scoreCol] = 0
    df = df.astype({scoreCol: "double"})
    df.sort_values(by = scoreCol, ascending = False, inplace = True)
    df.drop_duplicates(subset = [0,1,2], inplace = True)
    df.sort_index(inplace = True)
    os.makedirs(os.path.join(intDir, factor), exist_ok = True)
    df.to_csv(os.path.join(intDir, factor, bedname), sep = "\t", index = False, header = False)

with concurrent.futures.ProcessPoolExecutor(36) as executor:
    beds = glob.glob(normDir + "*/*.csv")
    list(tqdm(executor.map(intersect, beds), total = len(beds)))

def num(x):
    try:
        return int(x.split("/")[-1].split(".bed")[0].split("_")[0])
    except:
        return x

# For DE
""" factors = os.listdir(intDir)
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
cobind.to_csv(os.path.join(cobindDir, "DE.csv"), sep = "\t") """

# For ATAC
factors = os.listdir(intDir)
sortedBeds = []
for factor in factors:
    beds = glob.glob(os.path.join(intDir, factor, "*.csv"))
    beds.sort(key = num)
    sortedBeds.append(beds)

cobind = pd.DataFrame()
for i in range(len(sortedBeds)):
    factor = sortedBeds[i][0].split("/")[-2]
    temp = pd.read_csv(sortedBeds[i][0], sep = "\t")
    cobind[factor] = temp.iloc[:, -1].to_list()
os.makedirs(cobindDir, exist_ok = True)
cobind.to_csv(os.path.join(cobindDir, "atac.csv"), sep = "\t")