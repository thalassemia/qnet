import pandas as pd
import os

gtf = os.path.expandvars('$SCRATCH/data/gencode.v29.gtf')
gtf = pd.read_csv(gtf, comment = '#', header = None, sep = '\t')
transcripts = gtf.loc[gtf[2] == 'transcript']

info = [i.split('; ') for i in transcripts[8]]
transcript_id = [i[1].split('\"')[1] for i in info]
gene_name = [i[3].split('\"')[1] for i in info]
output = pd.DataFrame({'Transcript ID': transcript_id, 'HGNC symbol': gene_name})

tfs = os.path.expandvars("$SCRATCH/data/Database.csv")
tfs = pd.read_csv(tfs)
tfs = tfs.loc[tfs["Is TF?"] == "Yes"]

output = output.merge(tfs[["HGNC symbol"]], on = "HGNC symbol")

output.to_csv(os.path.expandvars('$SCRATCH/data/enst2hgnc.txt'), header = False, index = False, sep = '\t')