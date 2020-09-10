import pandas as pd
import os
from tqdm import tqdm
import concurrent.futures
from shutil import copy
import glob
import shlex
import subprocess
import io
import itertools
import json

# Adjustable settings
index_dir = os.path.expandvars('$SCRATCH/hmchip')
bed_dir = os.path.expandvars('$SCRATCH/hmchip/human_hm')
out_dir = os.path.expandvars('$SCRATCH/hmchip/h4k20me3')
gtf = os.path.expandvars('$SCRATCH/gencode.v35.annotation.gtf')
deTFtarget_dir = os.path.expandvars('$SCRATCH/output/detargets')
peakdistance = 0
cores = 8

# Get the IDs of all relevant bed files
os.makedirs(out_dir, exist_ok=True)
index = pd.read_csv(os.path.join(index_dir, 'human_hm_full_QC.txt'), sep='\t')
index = index.loc[index['Factor'] == 'H4K20me3'] 
index = index[index.columns.difference(['Species', 'GSMID', 'Factor'], sort=False)]

def get_file(id):
    file = glob.glob(os.path.join(bed_dir, str(id)) + '*')
    copy(os.path.join(bed_dir, ''.join(file)), out_dir)
    return file

def read_file(file):
    return pd.read_csv(file, sep = '\t', header = None)

with concurrent.futures.ProcessPoolExecutor(cores) as executor:
    # Isolate relevant bed files
    files = list(tqdm(executor.map(get_file, index.loc[:,'DCid'].to_numpy()), total=len(index.index)))
    files = list(itertools.chain.from_iterable(files))
    # Read from and concatenate all relevant bed files
    contents = list(tqdm(executor.map(read_file, files), total = len(files)))
    files = pd.concat(contents).sort_values([0,1,2])
    files.to_csv(os.path.join(out_dir, f'allpeaks{peakdistance}.csv'), sep = '\t', index = False, header = False)
    # Merge peaks within "peakdistance" base pairs into a single peak
    command = shlex.split(f"bedtools merge -d {peakdistance} -c 1 -o count -i {os.path.join(out_dir, f'allpeaks{peakdistance}.csv')}")
    result = subprocess.run(command, stdout=subprocess.PIPE)
    output = result.stdout
    df = pd.read_csv(io.BytesIO(output), sep="\t", header = None)
    df.to_csv(os.path.join(out_dir, f'uniquepeaks{peakdistance}.csv'), sep='\t', header = False, index = False)

# Use UROPA and Gencode v35 to annotate peaks with nearest protein coding gene
os.chdir(out_dir)
with open('query.json', 'w') as f:
    data = {}
    data['queries'] = []
    data['queries'].append({
        'feature':'gene',
        'feature.anchor':'start',
        'distance':3000,
        'filter.attribute':'gene_type',
        'attribute.value':'protein_coding'
    })
    data['show_attributes'] = 'gene_name'
    data['gtf'] = gtf
    data['bed'] = os.path.join(out_dir, f'uniquepeaks{peakdistance}.csv')
    data['threads'] = 8
    json.dump(data, f)
command = shlex.split(f"uropa -i query.json -s")
subprocess.run(command)

# Find target overlaps with DE TFs
peaks = pd.read_csv(f'uniquepeaks{peakdistance}_finalhits.txt', sep = '\t')
peaks = peaks.dropna()
peaks.to_csv(f'uniquepeaks{peakdistance}_finalhits.txt', sep = '\t')
def overlap(file):
    targets = pd.read_csv(os.path.join(deTFtarget_dir, file), sep = '\t')
    targets.rename(columns = {'GeneSymbol':'gene_name'}, inplace = True)
    merged = targets.merge(peaks, on = 'gene_name')
    merged = merged.iloc[:, [0,1,2,4,6,7,8,9,10,11,13,16,18,19]]
    merged['deTF'] = [str.split(file, ".")[0]] * len(merged.index)
    return merged
with concurrent.futures.ProcessPoolExecutor(cores) as executor:
    deTFtargets = os.listdir(deTFtarget_dir)
    merged = list(tqdm(executor.map(overlap, deTFtargets), total = len(deTFtargets)))
    merged = pd.concat(merged).sort_values(['deTF', 'peak_chr', 'peak_start', 'peak_end'])
    merged.to_csv(f'mergedtargets{peakdistance}.txt', sep = '\t', index = False)