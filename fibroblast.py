import concurrent.futures
import pandas as pd
from tqdm import tqdm
import os
import shutil
import glob
import subprocess
import shlex
import io

index_path = os.path.expandvars('$SCRATCH/tfchip/human_factor_full_QC.txt')
bed_dir = os.path.expandvars('$SCRATCH/tfchip/human_factor')
out_dir = os.path.expandvars('$SCRATCH/tfchip/skin_fibroblast')
peakdistance = 0
cores = 8
os.makedirs(out_dir, exist_ok = True)

index = pd.read_csv(index_path, sep = '\t')
index = index[index['Cell_type'] == 'Fibroblast']
index = index[index['Tissue_type'] == 'Skin']
index = index[index['FastQC'] >= 25]
index = index[index['UniquelyMappedRatio'] >= 0.5]
index = index[index['PBC'] >= 0.5]
index = index[index['PeaksUnionDHSRatio'] >= 0.7]
factors = index.loc[:,'Factor'].drop_duplicates()

def move(i):
    os.makedirs(os.path.join(out_dir, index.loc[i,'Factor']), exist_ok = True)
    bedfile = glob.glob(os.path.join(bed_dir, str(index.loc[i, 'DCid'])) + '*.bed')[0]
    dst = os.path.join(out_dir, index.loc[i, 'Factor'], str.split(bedfile,'/')[-1])
    shutil.copy(bedfile, dst)

def merge(factor):
    os.chdir(os.path.join(out_dir, factor))
    beds = []
    for bed in os.listdir():
        beds.append(pd.read_csv(bed, header = None, sep = '\t'))
    merged = pd.concat(beds)
    merged.iloc[:,[1,2]] = merged.iloc[:, [1,2]].astype('int64')
    merged = merged.sort_values([0,1,2])
    merged.to_csv(f'allpeaks{peakdistance}.csv', sep = '\t', header = False, index = False)
    command = shlex.split(f'bedtools merge -d {peakdistance} -c 1 -o count -i allpeaks{peakdistance}.csv')
    result = subprocess.run(command, stdout=subprocess.PIPE)
    output = result.stdout
    df = pd.read_csv(io.BytesIO(output), sep='\t', header = None)
    df.to_csv(f'uniquepeaks{peakdistance}.csv', sep='\t', header = False, index = False)

with concurrent.futures.ProcessPoolExecutor(cores) as executor:
    list(tqdm(executor.map(move, index.index), total = len(index.index)))
    list(tqdm(executor.map(merge, factors), total = len(factors)))