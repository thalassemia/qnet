import subprocess
import pandas as pd
import shlex
import io
import os
import glob
from pathlib import Path
import concurrent.futures
from tqdm import tqdm

peakDir = "/home/sean/peaks/"
outDir = "/home/sean/coassoc/"

def intersect(factor):
    cofactors = next(os.walk(peakDir))[1]
    for cofactor in cofactors:
        if cofactor == factor:
            continue
        for a in glob.glob(peakDir + factor + "/*.bed"):
            for b in glob.glob(peakDir + cofactor + "/*.bed"):
                command = shlex.split(f"/home/sean/bedtools intersect -a {a} -b {b} -c")
                result = subprocess.run(command, stdout=subprocess.PIPE)
                output = result.stdout.decode('utf-8')
                df = pd.read_csv(io.StringIO(output), sep="\t")
                fileIDs = (a.split('/')[-1]).split('_')[0] + "_" + cofactor + "_" + (b.split('/')[-1]).split('_')[0] + ".csv"
                p = Path(f'{outDir}/{factor}/')
                p.mkdir(exist_ok=True)
                df.to_csv(f"{p}/{fileIDs}", sep=",", index=False)

filesystem = os.walk(peakDir)
factors = next(filesystem)[1]

with concurrent.futures.ProcessPoolExecutor() as executor:
    tqdm(executor.map(intersect, factors), total=len(factors))