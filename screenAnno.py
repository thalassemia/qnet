import os
import pandas as pd
import glob
import re

cobindDir = "/u/scratch/s/seanchea/tfchip/normalizedAnnoSALL2/qval/"
os.chdir(cobindDir)

factors = os.listdir()

for factor in factors:
    dirFiles = glob.glob(factor + "/*.bed")
    dirFiles.sort(key=lambda f: int(re.sub('\D', '', f)))
    print(dirFiles)