import pandas as pd
import os
import glob
import concurrent.futures
from tqdm import tqdm

peakDir = "/u/scratch/s/seanchea/goodBeds/"
outDir = "/u/scratch/s/seanchea/normalized"
rankBy = 6 
if rankBy == 6:
    outDir += "s"
else:
    outDir += "q"

def norm(factor, peakDir):
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
    list(tqdm(executor.map(norm, factors, [peakDir]*len(factors)), total=len(factors)))