import pandas as pd
import os
import subprocess
import json
import shlex
import concurrent.futures
from tqdm import tqdm
import glob
import numpy as np

out_dir = os.path.expandvars("$SCRATCH/tfchip/uropa")
gtf = os.path.expandvars("$SCRATCH/data/gencode.v35.annotation.gtf")
bed_dir = os.path.expandvars("$SCRATCH/tfchip/goodBedsSALL2/")
normDir = os.path.expandvars("$SCRATCH/tfchip/Unormalized/")
intDir = os.path.expandvars("$SCRATCH/tfchip/Uintersect/")
cobindDir = os.path.expandvars("$SCRATCH/tfchip/Ucobind/")
deFile = os.path.expandvars("$SCRATCH/data/bsf.csv")
de = pd.read_csv(deFile, sep = "\t", header = None)
de = de.loc[de[2] < 0.05]
de.rename(columns = {0:"gene_name"}, inplace = True)
de.drop_duplicates(subset = "gene_name", inplace = True)

# Use UROPA and Gencode v35 to annotate peaks with nearest protein coding gene
def anno(bed):
    pieces = bed.split("/")
    out = os.path.join(out_dir, pieces[-2], pieces[-1])
    os.makedirs(out + "/", exist_ok = True)
    with open('query.json', 'w') as f:
        data = {}
        data['queries'] = []
        data['queries'].append({
            'feature':'gene',
            'feature.anchor':'start',
            'distance':3000,
            'filter_attribute':'gene_type',
            'attribute_value':'protein_coding'
        })
        data['show_attributes'] = 'gene_name'
        data['gtf'] = gtf
        data['bed'] = bed
        data['threads'] = 1
        json.dump(data, f)
    command = shlex.split(f"uropa -i query.json -s -o {out}")
    subprocess.run(command)

""" with concurrent.futures.ProcessPoolExecutor(36) as executor:
    beds = glob.glob(bed_dir + "*/*.bed")
    list(tqdm(executor.map(anno, beds), total = len(beds))) """

def norm(a):
    temp = pd.read_csv(a, sep = "\t")
    temp.sort_values(by="peak_score", ascending=False, inplace=True)
    temp.index = range(0,len(temp.index))
    temp["peak_score"] = (-1*(temp.index)+len(temp.index))/len(temp.index)
    temp.dropna(subset = ["gene_name"], inplace = True)
    bed = a.split("/")[-2]
    factor = a.split("/")[-3]
    os.makedirs(os.path.join(normDir, factor), exist_ok=True)
    temp.to_csv(os.path.join(normDir, factor, bed), sep="\t", index=False)

""" with concurrent.futures.ProcessPoolExecutor(36) as executor:
    beds = glob.glob(out_dir + "/*/*/*finalhits.txt")
    list(tqdm(executor.map(norm, beds), total=len(beds))) """

def intersect(a):
    bed = pd.read_csv(a, sep = "\t")
    merged = de.merge(bed, on = "gene_name", how = "left")
    merged = merged.loc[:, ["gene_name", "peak_score"]]
    merged.fillna(0, inplace = True)
    bedname = a.split("/")[-1]
    factor = a.split("/")[-2]
    os.makedirs(os.path.join(intDir, factor), exist_ok = True)
    merged.to_csv(os.path.join(intDir, factor, bedname), sep = "\t", index = False)

""" with concurrent.futures.ProcessPoolExecutor(36) as executor:
    beds = glob.glob(normDir + "*/*.bed")
    list(tqdm(executor.map(intersect, beds), total = len(beds))) """

factors = os.listdir(intDir)
sortedBeds = []
for factor in factors:
    beds = glob.glob(os.path.join(intDir, factor, "*.bed"))
    beds.sort(key = lambda x: int(x.split("/")[-1].split(".bed")[0].split("_")[0]))
    sortedBeds.append(beds)

cobindMat = pd.DataFrame()
deGenes = []
for i in range(len(sortedBeds)):
    factor = sortedBeds[i][0].split("/")[-2]
    temp = pd.read_csv(sortedBeds[i][0], sep = "\t")
    temp.drop_duplicates(subset = "gene_name", inplace = True, keep = "first", ignore_index = True)
    cobindMat[factor] = temp["peak_score"]
    deGenes = temp["gene_name"].to_list()
cobindMat.index = deGenes
os.makedirs(cobindDir, exist_ok = True)
cobindMat.to_csv(os.path.join(cobindDir, "train.csv"), sep = "\t")
print(cobindMat)

cobindMat = pd.DataFrame()
deGenes = []
for i in range(len(sortedBeds)):
    factor = sortedBeds[i][0].split("/")[-2]
    temp = pd.read_csv(sortedBeds[i][0], sep = "\t")
    if len(sortedBeds[i]) > 1:
        temp = pd.read_csv(sortedBeds[i][1], sep = "\t")
    temp.drop_duplicates(subset = "gene_name", inplace = True, keep = "first", ignore_index = True)
    cobindMat[factor] = temp["peak_score"]
    deGenes = temp["gene_name"].to_list()
cobindMat.index = deGenes
os.makedirs(cobindDir, exist_ok = True)
cobindMat.to_csv(os.path.join(cobindDir, "test.csv"), sep = "\t")
print(cobindMat)