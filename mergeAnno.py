import pandas as pd
import concurrent.futures
import os
from tqdm import tqdm

normBase = "/u/scratch/s/seanchea/tfchip/normalizedSALL2/qval"
os.chdir(normBase)
outBase = "/u/scratch/s/seanchea/tfchip/normalizedAnnoSALL2/qval"
annoBase = "/u/scratch/s/seanchea/tfchip/goodAnnoSALL2"

factors = os.listdir()

def merge(factor):
    normDir = os.path.join(normBase, factor)
    normBeds = os.listdir(normDir)
    normBeds.sort()
    annoDir = os.path.join(annoBase, factor)
    annoBeds = os.listdir(annoDir)
    annoBeds.sort()
    outDir = os.path.join(outBase, factor)
    os.makedirs(outDir, exist_ok=True)
    for i in range(len(normBeds)):
        ibed = pd.read_csv(os.path.join(normDir, normBeds[i]), sep = "\t", header = None)
        ibed.columns = ["chr", "TSS", "TTS", "V4", "Norm"]
        abed = pd.read_csv(os.path.join(annoDir, annoBeds[i]))
        merged = ibed.merge(abed, how="left", on=["V4"])
        merged = merged.loc[:, ["chr", "TSS", "TTS", "V4", "annotation", "geneId", "transcriptId", "distanceToTSS", "Norm"]]
        merged.to_csv(os.path.join(outDir, normBeds[i]), sep = "\t", index = False, header = False)


with concurrent.futures.ProcessPoolExecutor(36) as executor:
    list(tqdm(executor.map(merge, factors), total = len(factors)))