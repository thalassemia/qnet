import pandas as pd
import os
from shutil import copy
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr, fisher_exact

os.makedirs(os.path.expandvars("$SCRATCH/compareCond/plots/"), exist_ok = True)
tfPath = os.path.expandvars("$SCRATCH/data/Database.csv")
tf = pd.read_csv(tfPath)
tf.rename(columns = {"HGNC symbol": "gene_name"}, inplace = True)
tf['gene_name'] = tf['gene_name'].str.strip()
# merge using Ensembl IDs is more stable than with HGNC symbols
tf['Ensembl ID'] = tf['Ensembl ID'].str.strip()
tf = tf.loc[tf["Is TF?"] == "Yes"]

ciPath = os.path.expandvars("$SCRATCH/rna-seq/deseq_CIvsP_gencodev29_allgenes.csv")
ci = pd.read_csv(ciPath)
ci['gene_name'] = ci['gene_name'].str.strip()
# in all cases we want padj <= 0.05 and abs(log2FoldChange) >= 1
ci = ci.loc[ci["padj"] <= 0.05]
ci = ci.loc[abs(ci["log2FoldChange"]) >= 1]
ci["gene_id"] = [i.split('.')[0] for i in ci["gene_id"]]
ci.rename(columns = {"gene_id": "Ensembl ID"}, inplace = True)
ci_tfs = ci[["Ensembl ID", "gene_name", "log2FoldChange", "padj"]].merge(tf[["Ensembl ID"]], on = "Ensembl ID")
ci_tfs.to_csv(os.path.expandvars("$SCRATCH/compareCond/ci_tfs.csv"))

cirPath = os.path.expandvars("$SCRATCH/rna-seq/deseq_CIvsCIR_gencodev29_allgenes_122220.csv")
cir = pd.read_csv(cirPath)
cir['gene_name'] = cir['gene_name'].str.strip()
cir = cir.loc[cir["padj"] <= 0.05]
cir = cir.loc[abs(cir["log2FoldChange"]) >= 1]
cir["gene_id"] = [i.split('.')[0] for i in cir["gene_id"]]
cir.rename(columns = {"gene_id": "Ensembl ID"}, inplace = True)
cir_tfs = cir[["Ensembl ID", "gene_name", "log2FoldChange", "padj"]].merge(tf[["Ensembl ID"]], on = "Ensembl ID")
cir_tfs.to_csv(os.path.expandvars("$SCRATCH/compareCond/cir_tfs.csv"))

sirPath = os.path.expandvars("$SCRATCH/rna-seq/deseq_SSvsSSR_gencodev29_allgenes.csv")
sir = pd.read_csv(sirPath)
sir['gene_name'] = sir['gene_name'].str.strip()
sir = sir.loc[sir["padj"] <= 0.05]
sir = sir.loc[abs(sir["log2FoldChange"]) >= 1]
sir["gene_id"] = [i.split('.')[0] for i in sir["gene_id"]]
sir.rename(columns = {"gene_id": "Ensembl ID"}, inplace = True)
sir_tfs = sir[["Ensembl ID", "gene_name", "log2FoldChange", "padj"]].merge(tf[["Ensembl ID"]], on = "Ensembl ID")
sir_tfs.to_csv(os.path.expandvars("$SCRATCH/compareCond/sir_tfs.csv"))

ssPath = os.path.expandvars("$SCRATCH/data/deseq_SSvsP_gencodev29_allgenes_021120.txt")
ss = pd.read_csv(ssPath, sep = "\t")
ss['gene_name'] = ss['gene_name'].str.strip()
ss = ss.loc[ss["padj"] <= 0.05]
ss = ss.loc[abs(ss["log2FoldChange"]) >= 1]
ss["gene_id"] = [i.split('.')[0] for i in ss["gene_id"]]
ss.rename(columns = {"gene_id": "Ensembl ID"}, inplace = True)
ss_tfs = ss[["Ensembl ID", "gene_name", "log2FoldChange", "padj"]].merge(tf[["Ensembl ID"]], on = "Ensembl ID")
merged = ci[["gene_name", "log2FoldChange", "padj"]].merge(ss[["gene_name", "log2FoldChange", "padj"]], on = "gene_name", indicator = True, how = "outer")
merged.rename(columns = {"log2FoldChange_x": "l2FC_CIvP", "padj_x": "padj_CIvP", "log2FoldChange_y": "l2FC_SSvP", "padj_y": "padj_SSvP",}, inplace = True)
merged.to_csv(os.path.expandvars("$SCRATCH/compareCond/allDEgenes.csv"))

deGenes = [ss, ci, sir, cir]
outNames = ['ss', 'ci', 'sir', 'cir']
deTFs = [ss_tfs, ci_tfs, sir_tfs, cir_tfs]

gene = np.zeros((4,4))
tf = np.zeros((4,4))
fisher = np.zeros((4,4))
union = deTFs[0][["Ensembl ID", "gene_name", "log2FoldChange", "padj"]].merge(deTFs[0][["Ensembl ID", "gene_name", "log2FoldChange", "padj"]], on = "Ensembl ID", indicator = True, how = "outer")

for i in range(4):
    for j in range(i,4):
        merged = deGenes[i][["Ensembl ID", "gene_name", "log2FoldChange", "padj"]].merge(deGenes[j][["Ensembl ID", "gene_name", "log2FoldChange", "padj"]], on = "Ensembl ID", indicator = True, how = "outer")
        merged.to_csv(os.path.expandvars("$SCRATCH/compareCond/" + outNames[i] + "_" + outNames[j] + "_genes.csv"))
        # count the total number of DE genes in any of the DESeq2 outputs
        #union = union.merge(merged, on = "Ensembl ID", how = "outer")
        same = merged.loc[np.sign(merged.loc[:,'log2FoldChange_x']) == np.sign(merged.loc[:,'log2FoldChange_y'])]
        diff = merged.loc[np.sign(merged.loc[:,'log2FoldChange_x']) != np.sign(merged.loc[:,'log2FoldChange_y'])]
        # theoretically there should be no output here because i == j means we are comparing the same files
        if i == j:
            print(diff)
        # create plot showing correlation between the log 2 fold changes of DE genes
        plt.scatter(merged.loc[merged["_merge"] == "both"]['log2FoldChange_x'], merged.loc[merged["_merge"] == "both"]['log2FoldChange_y'], c = "steelblue", label = 'DE Genes', alpha = 0.5)
        coef = np.polyfit(merged.loc[merged["_merge"] == "both"]['log2FoldChange_x'], merged.loc[merged["_merge"] == "both"]['log2FoldChange_y'], 1)
        eq2 = np.poly1d(coef)
        r2, p2 = pearsonr(merged.loc[merged["_merge"] == "both"]['log2FoldChange_x'], merged.loc[merged["_merge"] == "both"]['log2FoldChange_y'])
        plt.plot(merged.loc[merged["_merge"] == "both"]['log2FoldChange_x'], eq2(merged.loc[merged["_merge"] == "both"]['log2FoldChange_x']), c = 'r', label = f'r = {round(r2, 3)}')
        # count the number of genes that appear in both compared lists and change in the same or different directions
        gene[j,i] = np.sum(diff['_merge'] == 'both')
        gene[i,j] = np.sum(same['_merge'] == 'both')
        # use fisher exact test to calculate p-value of observing the given overlap in DE genes (using 6432 total DE genes)
        if i != j:
            degCont = [[np.sum(merged['_merge'] == 'both'), np.sum(merged['_merge'] == 'left_only')], [np.sum(merged['_merge'] == 'right_only'), 6432 - np.sum(merged['_merge'] == 'both') - np.sum(merged['_merge'] == 'left_only') - np.sum(merged['_merge'] == 'right_only')]]
            odds, fisher[i,j] = fisher_exact(degCont, alternative = 'greater')
            print(degCont)
            print(fisher_exact(degCont, alternative = 'greater'))
        merged = deTFs[i][["Ensembl ID", "gene_name", "log2FoldChange", "padj"]].merge(deTFs[j][["Ensembl ID", "gene_name", "log2FoldChange", "padj"]], on = "Ensembl ID", indicator = True, how = "outer")
        merged.to_csv(os.path.expandvars("$SCRATCH/compareCond/" + outNames[i] + "_" + outNames[j] + "_tfs.csv"))
        # count the number of DE TFs in any of the DESeq2 outputs
        #union = union.merge(merged, on = "Ensembl ID", how = "outer")
        same = merged.loc[np.sign(merged.loc[:,'log2FoldChange_x']) == np.sign(merged.loc[:,'log2FoldChange_y'])]
        diff = merged.loc[np.sign(merged.loc[:,'log2FoldChange_x']) != np.sign(merged.loc[:,'log2FoldChange_y'])]
        # count the number of tfs that appear in both compared lists and change in the same or different directions
        tf[j,i] = np.sum(diff['_merge'] == 'both')
        tf[i,j] = np.sum(same['_merge'] == 'both')
        # use fisher exact test to calculate p-value of observing the given overlap in DE TFs (using 406 total DE TFs)
        if i != j:
            degCont = [[np.sum(merged['_merge'] == 'both'), np.sum(merged['_merge'] == 'left_only')], [np.sum(merged['_merge'] == 'right_only'), 406 - np.sum(merged['_merge'] == 'both') - np.sum(merged['_merge'] == 'left_only') - np.sum(merged['_merge'] == 'right_only')]]
            odds, fisher[j,i] = fisher_exact(degCont, alternative = 'greater')
            #print(degCont)
            #print(fisher_exact(degCont, alternative = 'greater'))
        # overlay plot showing correlation between the log 2 fold changes of DE TFs
        merged = merged.loc[merged['_merge'] == 'both']
        plt.scatter(merged["log2FoldChange_x"], merged["log2FoldChange_y"], c = 'darkorange', label = 'DE TFs', alpha = 0.5)
        plt.xlabel(outNames[i].upper() + " Log 2 Fold Change")
        plt.ylabel(outNames[j].upper() + " Log 2 Fold Change")
        coef = np.polyfit(merged["log2FoldChange_x"], merged['log2FoldChange_y'], 1)
        eq1 = np.poly1d(coef)
        r1, p1 = pearsonr(merged["log2FoldChange_x"], merged['log2FoldChange_y'])
        plt.plot(merged["log2FoldChange_x"], eq1(merged["log2FoldChange_x"]), c = 'darkorange', label = f'r = {round(r1, 3)}')
        plt.legend()
        plt.savefig(os.path.expandvars("$SCRATCH/compareCond/plots/" + outNames[i] + "_" + outNames[j] + ".png"))
        plt.clf()
np.savetxt(os.path.expandvars("$SCRATCH/compareCond/gene_counts.csv"), gene)
np.savetxt(os.path.expandvars("$SCRATCH/compareCond/tf_counts.csv"), tf)
np.savetxt(os.path.expandvars("$SCRATCH/compareCond/fisher.csv"), fisher)
#print(union.shape)