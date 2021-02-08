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

normDir = os.path.expandvars("$SCRATCH/tfchip/Cnormalized/")
intDir = os.path.expandvars("$SCRATCH/atac/Cintersect")
bedDir = os.path.expandvars("$SCRATCH/atac/Cbeds")
cobindDir = os.path.expandvars("$SCRATCH/atac/Ccobind/")
atac = pd.read_csv(os.path.expandvars("$SCRATCH/data/SS_P_diff_bound_sites.csv"))
atac = atac.loc[:, ["seqnames", "start", "end", "Fold", "FDR"]]
cleanPath = os.path.expandvars("$SCRATCH/atac/corr/atacAnno_clean.tsv")
atac.to_csv(cleanPath, sep = "\t", index = False, header = False)

def intersect(a):
    bedname = a.split("/")[-1]
    factor = a.split("/")[-2]
    bed = pd.read_csv(a, sep = "\t")
    bed = bed.loc[:, ["seqnames", "start", "end", "peak_score"]]
    os.makedirs(os.path.join(bedDir, factor), exist_ok = True)
    bed.to_csv(os.path.join(bedDir, factor, bedname), sep = "\t", index = False, header = False)
    command = shlex.split(f"bedtools intersect -a {cleanPath} -b {os.path.join(bedDir, factor, bedname)} -loj")
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

def cobind(dataName):
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
    cobind.to_csv(os.path.join(cobindDir, dataName + "_train.csv"), sep = "\t")

    cobind = pd.DataFrame()
    for i in range(len(sortedBeds)):
        factor = sortedBeds[i][0].split("/")[-2]
        temp = pd.read_csv(sortedBeds[i][0], sep = "\t")
        if len(sortedBeds[i]) > 1:
            temp = pd.read_csv(sortedBeds[i][1], sep = "\t")
        cobind[factor] = temp.iloc[:, -1].to_list()
    cobind.to_csv(os.path.join(cobindDir, dataName + "_test.csv"), sep = "\t")

cobind("chea")

normDir = os.path.expandvars("$SCRATCH/tfchip/Lnormalized/")
intDir = os.path.expandvars("$SCRATCH/atac/Lintersect/")
bedDir = os.path.expandvars("$SCRATCH/atac/Lbeds")
cobindDir = os.path.expandvars("$SCRATCH/atac/Lcobind/")

with concurrent.futures.ProcessPoolExecutor(36) as executor:
    beds = glob.glob(normDir + "*/*.csv")
    list(tqdm(executor.map(intersect, beds), total = len(beds)))

cobind("lisa")

normDir = os.path.expandvars("$SCRATCH/tfchip/Unormalized/")
intDir = os.path.expandvars("$SCRATCH/atac/Uintersect/")
bedDir = os.path.expandvars("$SCRATCH/atac/Ubeds")
cobindDir = os.path.expandvars("$SCRATCH/atac/Ucobind/")

with concurrent.futures.ProcessPoolExecutor(36) as executor:
    beds = glob.glob(normDir + "*/*.csv")
    list(tqdm(executor.map(intersect, beds), total = len(beds)))

cobind("up")