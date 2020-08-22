import pandas as pd
import os
import glob
import concurrent.futures
from tqdm import tqdm

peakDir = "/home/sean/tfchip/goodBeds/"
outDir = "/home/sean/tfchip/normalized"
rankBy = 6 # 8 for q-value (FDR-adjusted p-value), 6 for signal value (enrichment)
if(rankBy == 6):
    r = "q"
else:
    r = "s"

def norm(factor, peakDir):
    beds = glob.glob(peakDir + factor + "/*.bed")
    for a in beds:
        temp = pd.read_csv(a, sep="\t", header=None)
        temp = temp.iloc[:, [0,1,2,3,rankBy]]
        temp.sort_values(by=rankBy, ascending=False, inplace=True)
        temp.index = range(0,len(temp.index))
        temp[rankBy] = (-1*(temp.index)+len(temp.index))/len(temp.index)
        a = a.split("/")[-1]
        os.makedirs(os.path.join(outDir,factor,r), exist_ok=True)
        temp.to_csv(os.path.join(outDir,factor,r)+a, sep="\t", index=False, header=False)

factors = next(os.walk(peakDir))[1]
with concurrent.futures.ProcessPoolExecutor() as executor:
    list(tqdm(executor.map(norm, factors, [peakDir]*len(factors)), total=len(factors)))