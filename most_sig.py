import os
import concurrent.futures
from shutil import copy
import glob
import pandas as pd

indir = '/home/sean/BETA_output/'
outdir = '/home/sean/targets/'

def gene_name(file):
    filename = file.split('/')[4]
    gene = filename.split('_')[0]
    ud =  filename.split('_')[2]
    return (gene, ud)

def unique_targets(files, gene):
    targets = []
    for filename in files:
        df = pd.read_csv(filename, sep='\t', index_col=None, header=0)
        targets.append(df)
    frame = pd.concat(targets, axis=0, ignore_index=True)
    frame = frame.drop_duplicates('GeneSymbol')
    frame.sort_values(by='rank product', inplace=True)
    frame.to_csv(outdir + gene[0] + '_' + gene[1] + '.csv', sep='\t', index=False)

files = glob.glob(indir + '*.txt')

with concurrent.futures.ProcessPoolExecutor() as executor:
    genes = set(executor.map(gene_name, files))
    gene_files = []
    for gene in genes:
        gene_files.append(glob.glob(indir + gene[0] + '*' + gene[1]))
    executor.map(unique_targets, gene_files, genes)