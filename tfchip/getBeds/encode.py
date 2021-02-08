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

# compile dictionary of url's for highest quality peak files
beds = {}
good_outs = ['optimal IDR thresholded peaks', 'pseudoreplicated IDR thresholded peaks', 'conservative IDR thresholded peaks', 'IDR thresholded peaks']
for i in experiments:
    for j in i['files']:
        if j['file_format'] == 'bed' and j.get('assembly') == 'GRCh38' and j['output_type'] in good_outs:
                a = j['output_type'] == good_outs[3] and beds.get(j['dataset'],[0])[-1] not in good_outs[:3]
                b = j['output_type'] == good_outs[2] and beds.get(j['dataset'],[0])[-1] not in good_outs[:2]
                c = j['output_type'] == good_outs[1] and beds.get(j['dataset'],[0])[-1] not in good_outs[:1]
                d = j['output_type'] == good_outs[0]
                if a or b or c or d:
                    beds[j['dataset']] = [i['target']['label'], j['cloud_metadata']['url'], j['output_type']]

# method to enable multiprocessed downloading of files
def fetch(e):
    url = e[1]
    target = e[0]
    os.makedirs(downloads + target, exist_ok = True)
    os.system('cd ' + downloads + target + ' && wget {}'.format(url))

# turn dictionary into list for mapping
beds = [i for i in beds.values()]
with ProcessPoolExecutor(36) as executor:
    executor.map(fetch, beds)

# recursively extract files as needed
os.system('gunzip -r ' + downloads)