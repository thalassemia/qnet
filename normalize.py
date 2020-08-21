import pandas as pd
import os
import glob
import concurrent.futures
from tqdm import tqdm

peakDir = "/home/sean/tfchip/goodBeds/"

def norm(factor, peakDir):
    beds = glob.glob(peakDir + factor + "/*.bed")
    for a in beds:
        a_addIndex = pd.read_csv(a, sep="\t")
        a_addIndex['1'] = (-(a_addIndex.index+1)+len(a_addIndex.index))/len(a_addIndex.index)
        a_addIndex.to_csv(a, sep="\t", index=False)

factors = next(os.walk(peakDir))[1]
with concurrent.futures.ProcessPoolExecutor() as executor:
    list(tqdm(executor.map(norm, factors, [peakDir]*len(factors)), total=len(factors)))