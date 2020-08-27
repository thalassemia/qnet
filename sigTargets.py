import os
import concurrent.futures
from shutil import copy
import glob
import pandas as pd

indir = os.path.expandvars('$SCRATCH/output/betaMinus/')
outdir = os.path.expandvars('$SCRATCH/output/targets/')

def gene_name(file):
    filename = file.split('/')[-1]
    gene = filename.split('_')[0]
    # ud =  filename.split('_')[2]
    # return (gene, ud)
    return gene

def unique_targets(files, gene):
    targets = []
    for filename in files:
        df = pd.read_csv(filename, sep='\t', index_col=None, skip_blank_lines=False, skiprows=[0,1,2,3])
        targets.append(df)
    frame = pd.concat(targets, axis=0, ignore_index=True)
    frame = frame.drop_duplicates('GeneSymbol')
    frame.sort_values(by='Score', ascending=False, inplace=True)
    os.makedirs(outdir, exist_ok=True)
    #frame.to_csv(outdir + gene[0] + '_' + gene[1] + '.csv', sep='\t', index=False)
    frame.to_csv(outdir + gene + '.csv', sep='\t', index=False)

files = glob.glob(indir + '*targets.txt')

with concurrent.futures.ProcessPoolExecutor() as executor:
    genes = set(executor.map(gene_name, files))
    gene_files = []
    for gene in genes:
        #gene_files.append(glob.glob(indir + gene[0] + '*' + gene[1]))
        gene_files.append(glob.glob(indir + gene + '*targets.txt'))
    executor.map(unique_targets, gene_files, genes)