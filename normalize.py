import pandas as pd
import os
import glob
import concurrent.futures
from tqdm import tqdm

peakDir = os.path.expandvars("$SCRATCH/tfchip/goodBeds/")
outDir = os.path.expandvars("$SCRATCH/tfchip/normalized/")

def norm(factor):
    beds = glob.glob(peakDir + factor + "/*.bed")
    for a in beds:
        temp = pd.read_csv(a, sep="\t", header=None)
        temp = temp.iloc[:, [0,1,2,3,rankBy]]
        temp.sort_values(by=rankBy, ascending=False, inplace=True)
        temp.index = range(0,len(temp.index))
        temp[rankBy] = (-1*(temp.index)+len(temp.index))/len(temp.index)
        a = a.split("/")[-1]
        os.makedirs(os.path.join(outDir,factor), exist_ok=True)
        temp.to_csv(os.path.join(outDir,factor,a), sep="\t", index=False, header=False)

factors = next(os.walk(peakDir))[1]
with concurrent.futures.ProcessPoolExecutor() as executor:
    rankBy = 6
    sigOut = outDir + "signal"
    qOut = outDir + "q"
    outDir = sigOut
    list(tqdm(executor.map(norm, factors), total=len(factors)))
    rankBy = 8
    outDir = qOut
    list(tqdm(executor.map(norm, factors), total=len(factors)))