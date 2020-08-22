import subprocess
import pandas as pd
import shlex
import io
import os
import glob
from pathlib import Path
import concurrent.futures
from tqdm import tqdm

peakDir = "/home/sean/tfchip/normalized/"
outDir = "/home/sean/tfchip/intersectq/"

def intersect(factor):
    cofactors = next(os.walk(peakDir))[1]
    a = glob.glob(peakDir + factor + "/*.bed")[0]
    for cofactor in cofactors:
        b = glob.glob(peakDir + cofactor + "/*.bed")[0]
        command = shlex.split(f"/home/sean/bedtools intersect -a {a} -b {b} -loj")
        result = subprocess.run(command, stdout=subprocess.PIPE)
        output = result.stdout
        df = pd.read_csv(io.BytesIO(output), sep="\t", header = None)
        df = df.iloc[:,[0,1,2,3,9]]
        fileIDs = ((a.split('/')[-1]).split('.')[0]).split('*')[0] + "_" + cofactor + "_" + ((b.split('/')[-1]).split('.')[0]).split('*')[0] + ".csv"
        p = Path(f'{outDir}/{factor}/')
        os.makedirs(p, exist_ok=True)
        df.to_csv(f"{p}/{fileIDs}", index=False, header = False)
    return

filesystem = os.walk(peakDir)
factors = next(filesystem)[1]
with concurrent.futures.ProcessPoolExecutor() as executor:
    list(tqdm(executor.map(intersect, factors), total=len(factors)))