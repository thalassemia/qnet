import pandas as pd
import requests, json
import os
import numpy as np
from concurrent.futures import ProcessPoolExecutor

# Path to list of genes to query in ENCODE database (default: RNA-seq DE genes)
deTF = os.path.expandvars('$SCRATCH/data/deTF.csv')
deTF = pd.read_csv(deTF)

# directory to store downloaded files
downloads = os.path.expandvars('$SCRATCH/cobinding/beds/')
os.makedirs(downloads, exist_ok = True)

# fetch metadata for experiments that target deTFs and don't have read depth concerns
url = (
    'https://www.encodeproject.org/search/?type=Experiment'
    '&status=released'
    '&target.investigated_as=transcription+factor'
    '&replicates.library.biosample.treatments!=*'
    '&frame=embedded'
    '&format=json'
    '&assay_term_name=ChIP-seq'
    '&limit=all'
    '&assembly=GRCh38'
    '&files.output_type=IDR thresholded peaks'
    '&files.output_type=optimal IDR thresholded peaks'
    '&files.output_type=pseudoreplicated IDR thresholded peaks'
    '&audit.ERROR.category!=extremely low read depth'
)
for tf in deTF['gene_name']:
    url += '&target.label=' + tf

r = requests.get(url)
experiments = r.json()['@graph']

# directory to store downloaded files
downloads = os.path.expandvars('$SCRATCH/cobinding/beds/')
# compile dictionary of information for highest qualtiy peak files
metadata = {}

# convert information about type of file into number for easy comparison
def number(x):
    # prefer optimal IDR thresholded peaks over other two
    good_outs = ['conservative IDR thresholded peaks', 'pseudoreplicated IDR thresholded peaks', 'optimal IDR thresholded peaks']
    if x not in good_outs:
        return 0
    if x == good_outs[2]:
        return 3
    if x == good_outs[1]:
        return 2
    if x == good_outs[0]:
        return 1

for i in experiments:
    # get information about experiment target and cell type
    metadata[i['@id']] = [i['target']['label'], i["biosample_ontology"]["cell_slims"]] + [0]*11
    for j in i['files']:
        # ensure hg38 mapping
        if j.get('assembly') == 'GRCh38':
            if j.get('file_format') == 'bed':
                newType = number(j['output_type'])
                oldType = metadata[j['dataset']][-4]
                # number of replicates = length of replicate list
                newReps = len(j['biological_replicates'])
                oldReps = metadata[j['dataset']][-3]
                # only update metadata if more favorable type and/or more replicates included
                if (newType >= oldType) and (newReps >= oldReps):
                    # loop over all quality metric data to find FRiP
                    for k in j['quality_metrics']:
                        if 'frip' in k:
                            metadata[j['dataset']][6] = k.get('frip')
                    metadata[j['dataset']][9] = newType
                    metadata[j['dataset']][10] = newReps
                    # store downloaded path to make referencing file info easier
                    metadata[j['dataset']][11] = os.path.join(downloads, i['target']['label'], j['accession'] + '.bed')
                    # url to download from
                    metadata[j['dataset']][12] = j['cloud_metadata']['url']
            # look through all quality metrics for alignment files to find NRF, PBC1, PBC2, and RSC
            if j['output_type'] == 'alignments':
                # get read length
                metadata[j['dataset']][7] = j["mapped_read_length"]
                for k in j['quality_metrics']:
                    if 'NRF' in k:
                        metadata[j['dataset']][2] = k['NRF']
                    if 'PBC1' in k:
                        metadata[j['dataset']][3] = k['PBC1']
                    if 'PBC2' in k:
                        metadata[j['dataset']][4] = k['PBC2']
                    if 'RSC' in k:
                        metadata[j['dataset']][5] = k['RSC']
                    # get number of fragments
                    if 'distinct_fragments' in k:
                        metadata[j['dataset']][8] = k['distinct_fragments']
                    # if # of fragments not available, use # of filtered reads
                    elif 'processing_stage' in k and metadata[j['dataset']][8] == 0:
                        if k['processing_stage'] == 'filtered':
                            try:
                                metadata[j['dataset']][8] = k['total_reads']
                            except:
                                metadata[j['dataset']][8] = k['total']

metadataDF = pd.DataFrame.from_dict(metadata, orient = 'index')
metadataDF.rename(columns = {0: 'Target', 1: 'Cell Type(s)', 2: 'NRF', 3: 'PBC1', 4: 'PBC2', 5: 'RSC', 6: 'FRiP', 7: 'Read Length', 8: 'Total Usable Fragments', 9:'File Type', 10: 'Replicates', 11: 'Path', 12: 'URL'}, inplace = True)
# only keep highest quality files for each factor that meet ENCODE's most stringent quality criteria
metadataDF = metadataDF.loc[(metadataDF['NRF'] >= 0.9) & (metadataDF['PBC1'] >= 0.9) & (metadataDF['PBC2'] >= 10) & (metadataDF['Read Length'] >= 50) & (metadataDF['Total Usable Fragments'] >= 20000000), ].sort_values(by = 'FRiP', ascending = False)
metadataDF.to_csv('/u/scratch/s/seanchea/cobinding/beds/encodeKey.csv')
# method to enable multiprocessed downloading of files
def fetch(e):
    url = e.loc['URL']
    target = e.loc['Target']
    print(url, target)
    os.makedirs(downloads + target, exist_ok = True)
    os.system('cd ' + downloads + target + ' && wget {}'.format(url))

# turn dictionary into list for mapping
beds = [metadataDF.loc[i,] for i in metadataDF.index]
with ProcessPoolExecutor(36) as executor:
    executor.map(fetch, beds)

# recursively extract files as needed
os.system('gunzip -r ' + downloads)