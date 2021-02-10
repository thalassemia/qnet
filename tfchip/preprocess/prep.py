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

anno_dir = os.path.expandvars("$SCRATCH/cobinding/annotated/")
normDir = os.path.expandvars("$SCRATCH/cobinding/normalized/")
intDir = os.path.expandvars("$SCRATCH/cobinding/intersect/")
cobindDir = os.path.expandvars("$SCRATCH/cobinding/cobind/")
deFile = os.path.expandvars("$SCRATCH/data/deseq_SSvsP_gencodev29_allgenes_021120.txt")
de = pd.read_csv(deFile, sep = "\t")
de.rename(columns = {"gene_id": "Ensembl ID"}, inplace = True)
de = de.loc[:,['Ensembl ID', 'gene_name','log2FoldChange','padj']]
de['gene_name'] = de['gene_name'].str.strip()
# remove anything whitespace and version number (decimal) from Ensembl ID to keep things intercompatible
de['Ensembl ID'] = de['Ensembl ID'].str.strip()
de['Ensembl ID'] = [i.split('.')[0] for i in de['Ensembl ID']]
de = de.loc[de['padj'] <= 0.05]

# get separate lists of significant upregulated and downregulated genes
up = de.loc[de['log2FoldChange'] >= 0]
down = de.loc[de['log2FoldChange'] <= 0]

def norm(a):
    temp = pd.read_csv(a)
    # V9 is q-value (for Cistrome peaks) and V7 is enrichment signal (for ENCODE peaks)
    temp.sort_values(by=["V9", "V7"], ascending=True, inplace=True)
    # maximum possible score is 1 and lowest possible score is just greater than 0 (prevent erroneous clustering)
    temp.index = range(1,len(temp.index)+1)
    temp["peak_score"] = temp.index/len(temp.index)
    bed = a.split("/")[-1]
    factor = a.split("/")[-2]
    os.makedirs(os.path.join(normDir, factor), exist_ok=True)
    temp.to_csv(os.path.join(normDir, factor, bed), sep="\t", index=False)

""" with concurrent.futures.ProcessPoolExecutor(36) as executor:
    beds = glob.glob(anno_dir + "*/*.bed")
    list(tqdm(executor.map(norm, beds), total=len(beds))) """

def intersect(a):
    bed = pd.read_csv(a, sep = "\t")
    bed.rename(columns = {"geneId": "Ensembl ID"}, inplace = True)
    bed['Ensembl ID'] = [i.split('.')[0] for i in bed['Ensembl ID']]
    # find peaks annotated to promoters of DE genes
    bed = bed.loc[[i.find("Promoter")!=-1 for i in bed["annotation"]],]
    merged = de.merge(bed, on = "Ensembl ID", how = "left")
    merged = merged.loc[:, ["gene_name", "peak_score", "Ensembl ID"]]
    # only keep the highest scoring peak for each Ensembl ID
    merged.sort_values(by = "peak_score", ascending = False, inplace = True)
    merged.drop_duplicates(subset = "Ensembl ID", inplace = True)
    merged.sort_values(by = "Ensembl ID", inplace = True)
    # all non-matched DE promoters get auto-assigned a score of 0
    merged.fillna(0, inplace = True)
    bedname = a.split("/")[-1]
    factor = a.split("/")[-2]
    os.makedirs(os.path.join(intDir, factor), exist_ok = True)
    merged.to_csv(os.path.join(intDir, factor, bedname), sep = "\t", index = False)

with concurrent.futures.ProcessPoolExecutor(36) as executor:
    beds = glob.glob(normDir + "*/*.bed")
    list(tqdm(executor.map(intersect, beds), total = len(beds)))

# get rank of bed file (for Cistrome) or assign rank of 0 (preference to ENCODE files)
def num(x):
    try:
        return int(x.split("/")[-1].split(".bed")[0])
    except:
        return 0

factors = os.listdir(intDir)
sortedBeds = []
for factor in factors:
    beds = glob.glob(os.path.join(intDir, factor, "*.bed"))
    beds.sort(key = num)
    sortedBeds.append(beds)

# compile cobinding maps for different conditions
cobindUU = pd.DataFrame() # upregulated TFs at upregulated DE genes
cobindUD = pd.DataFrame() # upregulated TFs at downregulated DE genes
cobindDD = pd.DataFrame() # downregulated TFs at downregulated DE genes
cobindDU = pd.DataFrame() # downregulated TFs at upregulated DE genes
cobindUp = pd.DataFrame() # upregulated TFs at all DE genes
cobindDown = pd.DataFrame() # downregulated TFs at all DE genes
cobind_Up = pd.DataFrame() # all DE TFs at upregulated DE genes
cobind_Down = pd.DataFrame() # all DE TFs at downregulated DE genes
cobindAll = pd.DataFrame() # all DE TFs at all DE genes
deGenes = []

for i in sortedBeds:
    factor = i[0].split("/")[-2]
    temp = pd.read_csv(i[0], sep = "\t")
    temp.sort_values(by = "Ensembl ID", inplace = True)
    up_genes = temp["Ensembl ID"].isin(up["Ensembl ID"])
    if factor in up["gene_name"].to_list():
        cobindUU[factor] = temp.loc[up_genes, "peak_score"]
        cobindUD[factor] = temp.loc[~up_genes, "peak_score"]
        cobindUp[factor] = temp["peak_score"]
    else:
        cobindDU[factor] = temp.loc[up_genes, "peak_score"]
        cobindDD[factor] = temp.loc[~up_genes, "peak_score"]
        cobindDown[factor] = temp["peak_score"]
    upDEG = temp.loc[up_genes, "gene_name"].to_list()
    downDEG = temp.loc[~up_genes, "gene_name"].to_list()
    allDEG = temp["gene_name"].to_list()
    cobind_Up[factor] = temp.loc[up_genes, "peak_score"]
    cobind_Down[factor] = temp.loc[~up_genes, "peak_score"]
    cobindAll[factor] = temp["peak_score"]

cobindUU.index = upDEG
cobindDU.index = upDEG
cobindUD.index = downDEG
cobindDD.index = downDEG
cobindUp.index = allDEG
cobindDown.index = allDEG
cobind_Up.index = upDEG
cobind_Down.index = downDEG
cobindAll.index = allDEG

os.makedirs(cobindDir, exist_ok = True)
cobindUU.to_csv(os.path.join(cobindDir, "UU_train.csv"), sep = "\t")
cobindUD.to_csv(os.path.join(cobindDir, "DD_train.csv"), sep = "\t")
cobindDD.to_csv(os.path.join(cobindDir, "DD_train.csv"), sep = "\t")
cobindDU.to_csv(os.path.join(cobindDir, "DD_train.csv"), sep = "\t")
cobindUp.to_csv(os.path.join(cobindDir, "up_train.csv"), sep = "\t")
cobindDown.to_csv(os.path.join(cobindDir, "down_train.csv"), sep = "\t")
cobind_Up.to_csv(os.path.join(cobindDir, "_up_train.csv"), sep = "\t")
cobind_Down.to_csv(os.path.join(cobindDir, "_down_train.csv"), sep = "\t")
cobindAll.to_csv(os.path.join(cobindDir, "all_train.csv"), sep = "\t")


# compile cobinding maps for different conditions
cobindUU = pd.DataFrame() # upregulated TFs at upregulated DE genes
cobindUD = pd.DataFrame() # upregulated TFs at downregulated DE genes
cobindDD = pd.DataFrame() # downregulated TFs at downregulated DE genes
cobindDU = pd.DataFrame() # downregulated TFs at upregulated DE genes
cobindUp = pd.DataFrame() # upregulated TFs at all DE genes
cobindDown = pd.DataFrame() # downregulated TFs at all DE genes
cobind_Up = pd.DataFrame() # all DE TFs at upregulated DE genes
cobind_Down = pd.DataFrame() # all DE TFs at downregulated DE genes
cobindAll = pd.DataFrame() # all DE TFs at all DE genes
deGenes = []

for i in sortedBeds:
    factor = i[0].split("/")[-2]
    temp = pd.read_csv(i[0], sep = "\t")
    # use second highest ranked file if available to create "test" set
    if len(i) > 1:
        temp = pd.read_csv(i[1], sep = "\t")
    temp.sort_values(by = "Ensembl ID", inplace = True)
    up_genes = temp["Ensembl ID"].isin(up["Ensembl ID"])
    if factor in up["gene_name"].to_list():
        cobindUU[factor] = temp.loc[up_genes, "peak_score"]
        cobindUD[factor] = temp.loc[[not i for i in up_genes], "peak_score"]
        cobindUp[factor] = temp["peak_score"]
    else:
        cobindDU[factor] = temp.loc[up_genes, "peak_score"]
        cobindDD[factor] = temp.loc[[not i for i in up_genes], "peak_score"]
        cobindDown[factor] = temp["peak_score"]
    upDEG = temp.loc[up_genes, "gene_name"].to_list()
    downDEG = temp.loc[[not i for i in up_genes], "gene_name"].to_list()
    allDEG = temp["gene_name"].to_list()
    cobind_Up[factor] = temp.loc[up_genes, "peak_score"]
    cobind_Down[factor] = temp.loc[[not i for i in up_genes], "peak_score"]
    cobindAll[factor] = temp["peak_score"]

cobindUU.index = upDEG
cobindDU.index = upDEG
cobindUD.index = downDEG
cobindDD.index = downDEG
cobindUp.index = allDEG
cobindDown.index = allDEG
cobind_Up.index = upDEG
cobind_Down.index = downDEG
cobindAll.index = allDEG

os.makedirs(cobindDir, exist_ok = True)
cobindUU.to_csv(os.path.join(cobindDir, "UU_test.csv"), sep = "\t")
cobindUD.to_csv(os.path.join(cobindDir, "DD_test.csv"), sep = "\t")
cobindDD.to_csv(os.path.join(cobindDir, "DD_test.csv"), sep = "\t")
cobindDU.to_csv(os.path.join(cobindDir, "DD_test.csv"), sep = "\t")
cobindUp.to_csv(os.path.join(cobindDir, "up_test.csv"), sep = "\t")
cobindDown.to_csv(os.path.join(cobindDir, "down_test.csv"), sep = "\t")
cobind_Up.to_csv(os.path.join(cobindDir, "_up_test.csv"), sep = "\t")
cobind_Down.to_csv(os.path.join(cobindDir, "_down_test.csv"), sep = "\t")
cobindAll.to_csv(os.path.join(cobindDir, "all_test.csv"), sep = "\t")

# output a list of DE gene promoters that had no annotated TF binding
noBind = cobindAll[(cobindAll==0).all(axis = 1)]
with open(os.path.join(cobindDir, "noBinding.csv"), "w") as out:
    out.writelines([i + '\n' for i in noBind.index])