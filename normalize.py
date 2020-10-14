import pandas as pd
import os
import glob
import concurrent.futures
from tqdm import tqdm

peakDir = os.path.expandvars("$SCRATCH/tfchip/goodBedsSALL2/")
outDir = os.path.expandvars("$SCRATCH/tfchip/normalizedSALL2/")
cores = 36

def norm(factor):
    beds = glob.glob(peakDir + factor + "/*.bed")
    for a in beds:
        temp = pd.read_csv(a, sep = "\t", header = None)
        temp = temp.iloc[:, [0,1,2,3,rankBy]]
        temp.sort_values(by=rankBy, ascending=False, inplace=True)
        temp.index = range(0,len(temp.index))
        temp[rankBy] = (-1*(temp.index)+len(temp.index))/len(temp.index)
        a = a.split("/")[-1]
        os.makedirs(os.path.join(outDir,factor), exist_ok=True)
        temp.to_csv(os.path.join(outDir,factor,a), sep="\t", index=False, header=False)

factors = next(os.walk(peakDir))[1]
print(factors)
rankBy = 6
sigOut = outDir + "signal"
qOut = outDir + "qval"
outDir = sigOut
with concurrent.futures.ProcessPoolExecutor(cores) as executor:
    list(tqdm(executor.map(norm, factors), total=len(factors)))
rankBy = 8
outDir = qOut
with concurrent.futures.ProcessPoolExecutor(cores) as executor:
    list(tqdm(executor.map(norm, factors), total=len(factors)))