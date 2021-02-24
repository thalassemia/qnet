import os
import pandas as pd
import numpy as np
from lisa import FromGenes, FromRegions
from scipy.stats import fisher_exact
import math

degs = pd.read_csv('/u/scratch/s/seanchea/data/deseq_SSvsP_gencodev29_allgenes_021120.txt', sep = '\t')
degs = degs.loc[degs['padj'] <= 0.05,]
degs['gene_name'].to_csv('/u/scratch/s/seanchea/data/DEGs.txt', header = False, index = False)
degs = degs['gene_name'].tolist()

results, metadata = FromRegions.using_macs_output('hg38', degs, '/u/scratch/s/seanchea/data/ATACraw/SS/ss123_peaks.xls')
results_df = pd.DataFrame(results.to_dict())
results_df.to_csv('/u/scratch/s/seanchea/lisaFR.csv', index = False)
genes = FromGenes('hg38')
resultsG, metadataG = genes.predict(degs, num_background_genes=10000)
resultsG_df = pd.DataFrame(resultsG.to_dict())
resultsG_df.to_csv('/u/scratch/s/seanchea/lisaFG.csv', index = False)

results_df = pd.read_csv('/u/scratch/s/seanchea/lisaFR.csv')
resultsG_df = pd.read_csv('/u/scratch/s/seanchea/lisaFG.csv')

results_dedup = results_df.drop_duplicates(subset = "factor")
results_dedup = results_dedup.loc[results_dedup["summary_p_value"] <= 0.05,]
results_dedup.loc[:,"Rank"] = range(1, results_dedup.shape[0]+1)

resultsG_dedup = resultsG_df.drop_duplicates(subset = "factor")
resultsG_dedup = resultsG_dedup.loc[resultsG_dedup["summary_p_value"] <= 0.05,]
resultsG_dedup.loc[:,"Rank"] = range(1, resultsG_dedup.shape[0]+1)


merged = resultsG_dedup.merge(results_dedup, on = "factor", how = "outer", indicator = True)

cutoffs = [10, 50, 250, 500]

for i in cutoffs:
    both = merged.loc[(merged["Rank_x"] <= i) & (merged["Rank_y"] <= i),].shape[0]
    G_only = merged.loc[(merged["Rank_x"] <= i) & (merged["_merge"] == "left_only"),].shape[0]
    R_only = merged.loc[(merged["Rank_y"] <= i) & (merged["_merge"] == "right_only"),].shape[0]
    print(both, G_only, R_only)
    print(-math.log10(fisher_exact([[both, R_only], [G_only, 963-both-2*R_only]], alternative="greater")[1]))

chea = pd.read_csv('/u/scratch/s/seanchea/chea.tsv', sep = '\t')
chea.rename(columns={"TF": "factor"}, inplace = True)
#lisa_chea = results_dedup.merge(chea, on = "factor", how = "outer", indicator = True)
# to see number of genes shared between all lists
lisa_chea = merged.loc[:,["Rank_x", "Rank_y", "factor"]].merge(chea, on = "factor", how = "outer", indicator = True)
for i in cutoffs:
    #both = lisa_chea.loc[(lisa_chea["Rank_x"] <= i) & (lisa_chea["Rank_y"] <= i),].shape[0]
    both = lisa_chea.loc[(lisa_chea["Rank_x"] <= i) & (lisa_chea["Rank_y"] <= i) & (lisa_chea["Rank"] <= i),].shape[0]
    l_only = i- both
    c_only = i - both
    print(both, l_only, c_only)
    # Use 1964 for FromGenes, 1972 for FromRegions, 1980 for both
    print(-math.log10(fisher_exact([[both, l_only], [c_only, 1980-both-2*l_only]], alternative="greater")[1]))

detfs = pd.read_csv('/u/scratch/s/seanchea/data/deTF.csv')
detfs.rename(columns={"gene_name": "factor"}, inplace = True)
detf_chea = detfs.merge(chea, on = "factor", how = "outer", indicator = True)
for i in cutoffs:
    both = detf_chea.dropna()
    both = both.loc[both["Rank"] <= i, ].shape[0]
    c_only = i - both
    d_only = 220 - both
    print(both, c_only, d_only)
    print(-math.log10(fisher_exact([[both, c_only], [d_only, 1634-both-l_only-d_only]], alternative="greater")[1]))

detf_lisa = detfs.merge(resultsG_dedup, on = "factor", how = "outer", indicator = True)
for i in cutoffs:
    both = detf_lisa.dropna()
    both = both.loc[both["Rank"] <= i, ].shape[0]
    l_only = i - both
    d_only = 220 - both
    print(both, l_only, d_only)
    # Use 1037 for FromGenes, 1074 for FromRegions
    print(-math.log10(fisher_exact([[both, l_only], [d_only, 1037-both-l_only-d_only]], alternative="greater")[1]))