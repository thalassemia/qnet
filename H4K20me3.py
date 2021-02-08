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
gtf = os.path.expandvars('$SCRATCH/data/gencode.v29.gtf')
deTFtarget_dir = os.path.expandvars('$SCRATCH/outputNoSALL2/detargets')
deData = os.path.expandvars('$SCRATCH/data/deseq_SSvsP_gencodev29_allgenes_021120.txt')
peakdistance = 0
cores = 36

# Get the IDs of all relevant bed files
os.makedirs(out_dir, exist_ok=True)
index = pd.read_csv(os.path.join(index_dir, 'human_hm_full_QC.txt'), sep='\t')
index = index.loc[index['Factor'] == 'H4K20me3'] 
index = index[index.columns.difference(['Species', 'GSMID', 'Factor'], sort=False)]

def get_file(id):
    fileName = glob.glob(os.path.join(bed_dir, str(id)) + '_*')[0]
    copy(os.path.join(bed_dir, fileName), os.path.join(out_dir, 'beds'))
    return fileName

def read_file(fileName):
    return pd.read_csv(fileName, sep = '\t', header = None)

with concurrent.futures.ProcessPoolExecutor(cores) as executor:
    # Isolate relevant bed files
    files = list(tqdm(executor.map(get_file, index.loc[:,'DCid'].to_numpy()), total=len(index.index)))
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

# Use UROPA and GENCODE V29 to annotate peaks with nearest protein coding gene
os.chdir(out_dir)
with open('query.json', 'w') as f:
    data = {}
    data['queries'] = []
    data['queries'].append({
        'feature':'gene',
        'feature.anchor':'start',
        'distance':3000,
        'filter_attribute':'gene_type',
        'attribute_value':'protein_coding'
    })
    data['show_attributes'] = 'gene_name'
    data['gtf'] = gtf
    data['bed'] = os.path.join(out_dir, f'uniquepeaks{peakdistance}.csv')
    data['threads'] = cores
    json.dump(data, f)
command = shlex.split(f"uropa -i query.json -s")
subprocess.run(command)

# Find target overlaps with DE TFs
peaks = pd.read_csv(f'uniquepeaks{peakdistance}_finalhits.txt', sep = '\t')
peaks = peaks.dropna()
peaks.to_csv(f'uniquepeaks{peakdistance}_finalhits.txt', sep = '\t', index = False)
def overlap(file):
    targets = pd.read_csv(os.path.join(deTFtarget_dir, file), sep = ',')
    targets.rename(columns = {'GeneSymbol':'gene_name'}, inplace = True)
    merged = targets.merge(peaks, on = 'gene_name')
    merged.rename(columns = {'#Chromsome':'Chromosome'}, inplace = True)
    merged = merged.loc[:, ['Chromosome', 'TSS', 'TTS', 'Score', 'Strand', 'gene_name',
        'log2FoldChange', 'padj', 'peak_chr', 'peak_start', 'peak_end',
        'distance', 'relative_location', 'feat_ovl_peak', 'peak_ovl_feat']]
    merged['deTF'] = [str.split(file, ".")[0]] * len(merged.index)
    return merged
with concurrent.futures.ProcessPoolExecutor(cores) as executor:
    deTFtargets = os.listdir(deTFtarget_dir)
    merged = list(tqdm(executor.map(overlap, deTFtargets), total = len(deTFtargets)))
    merged = pd.concat(merged, join='inner').sort_values(['deTF', 'peak_chr', 'peak_start', 'peak_end'])
    merged.to_csv(f'mergedtargets{peakdistance}.txt', sep = '\t', index = False)

# Merge H4K20me3-annotated genes with quiescence DE data to find patterns
de = pd.read_csv(deData, sep = '\t')
de = de.loc[de.padj <= 0.05]
merged = de.merge(peaks, on = 'gene_name')
merged.rename(columns = {'#Chromsome':'Chromosome'}, inplace = True)
merged = merged.loc[:, ['peak_chr', 'peak_start', 'peak_end', 'feat_start', 'feat_end', 'feat_strand',
    'distance', 'relative_location', 'feat_ovl_peak', 'peak_ovl_feat', 'gene_name', 'log2FoldChange', 'padj']]
merged.to_csv(f'allDE{peakdistance}.txt', sep = '\t', index = False)