import subprocess
import pandas as pd
import shlex
import io
import os
import glob
from pathlib import Path
import concurrent.futures
from tqdm import tqdm

peakDir = os.path.expandvars('$SCRATCH/tfchip/normalized/')
outDir = os.path.expandvars('$SCRATCH/output/intersect/')

def intersect(factor):
    cofactors = next(os.walk(peakDir))[1]
    a = glob.glob(peakDir + factor + "/*.bed")[0]
    for cofactor in cofactors:
        b = glob.glob(peakDir + cofactor + "/*.bed")[0]
        command = shlex.split(f"bedtools intersect -a {a} -b {b} -loj")
        result = subprocess.run(command, stdout=subprocess.PIPE)
        output = result.stdout
        df = pd.read_csv(io.BytesIO(output), sep="\t", header = None)
        df = df.iloc[:,[0,1,2,3,9]]
        fileIDs = ((a.split('/')[-1]).split('.')[0]).split('*')[0] + "_" + cofactor + "_" + ((b.split('/')[-1]).split('.')[0]).split('*')[0] + ".csv"
        p = Path(f'{outDir}/{factor}/')
        os.makedirs(p, exist_ok=True)
        df.to_csv(f"{p}/{fileIDs}", index=False, header = False)
    return

sigIn = peakDir + "signal/"
sigOut = outDir + "signal/"
qIn = peakDir + "qval/"
qOut = outDir + "qval/"
peakDir = sigIn
outDir = sigOut
with concurrent.futures.ProcessPoolExecutor() as executor:
    factors = next(os.walk(peakDir))[1]
    list(tqdm(executor.map(intersect, factors), total=len(factors)))
peakDir = qIn
outDir = qOut
with concurrent.futures.ProcessPoolExecutor() as executor:
    factors = next(os.walk(peakDir))[1]
    list(tqdm(executor.map(intersect, factors), total=len(factors)))