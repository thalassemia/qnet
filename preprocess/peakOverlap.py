import subprocess
import pandas as pd
import shlex
import io
import os
import glob
from pathlib import Path
import concurrent.futures
from tqdm import tqdm
import re

peakDir = os.path.expandvars('$SCRATCH/tfchip/normalizedAnnoSALL2/')
outDir = os.path.expandvars('$SCRATCH/tfchip/intersectAnnoSALL2/')
cores = 36

def intersect(factor):
    cofactors = next(os.walk(peakDir))[1]
    files = glob.glob(peakDir + factor + "/*.bed")
    files.sort(key=lambda f: int(re.sub('\D', '', f)))
    for i in range(2):
        if len(files) > 1:
            a = files[i]
        else:
            a = files[0]
        for cofactor in cofactors:
            cofiles = glob.glob(peakDir + cofactor + "/*.bed")
            cofiles.sort(key=lambda f: int(re.sub('\D', '', f)))
            if len(cofiles) > 1:
                b = cofiles[i]
            else:
                b = cofiles[0]
            command = shlex.split(f"bedtools intersect -a {a} -b {b} -loj")
            result = subprocess.run(command, stdout=subprocess.PIPE)
            output = result.stdout
            df = pd.read_csv(io.BytesIO(output), sep="\t", header = None)
            df = df.iloc[:,[0,1,2,3,4,5,6,7,-1]]
            df.columns = ["chr", "start", "end", "peak", "annnotation", "geneId", "transcriptId", "distanceToTSS", "norm_score"]
            fileIDs = ((a.split('/')[-1]).split('.')[0]) + "_" + cofactor + "_" + ((b.split('/')[-1]).split('.')[0]) + ".csv"
            p = Path(f'{outDir}/{factor}/')
            os.makedirs(p, exist_ok=True)
            df.to_csv(f"{p}/{fileIDs}", index=False, header = True)
    return

sigIn = peakDir + "signal/"
sigOut = outDir + "signal/"
qIn = peakDir + "qval/"
qOut = outDir + "qval/"
""" peakDir = sigIn
outDir = sigOut
with concurrent.futures.ProcessPoolExecutor(cores) as executor:
    factors = next(os.walk(peakDir))[1]
    list(tqdm(executor.map(intersect, factors), total=len(factors))) """
peakDir = qIn
outDir = qOut
with concurrent.futures.ProcessPoolExecutor(cores) as executor:
    factors = next(os.walk(peakDir))[1]
    list(tqdm(executor.map(intersect, factors), total=len(factors)))