import pandas as pd
import os
from shutil import copy
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import pearsonr, hypergeom

os.makedirs(os.path.expandvars("$SCRATCH/compareCond/"), exist_ok = True)
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

same_gene = np.zeros((4,4))
diff_gene = np.zeros((4,4))
same_tf = np.zeros((4,4))
diff_tf = np.zeros((4,4))

for i in range(4):
    for j in range(i,4):
        merged = deGenes[i][["Ensembl ID", "gene_name", "log2FoldChange", "padj"]].merge(deGenes[j][["Ensembl ID", "gene_name", "log2FoldChange", "padj"]], on = "Ensembl ID", indicator = True, how = "outer")
        merged.to_csv(os.path.expandvars("$SCRATCH/compareCond/" + outNames[i] + "_" + outNames[j] + "_genes.csv"))
        same = merged.loc[np.sign(merged.loc[:,'log2FoldChange_x']) == np.sign(merged.loc[:,'log2FoldChange_y'])]
        diff = merged.loc[np.sign(merged.loc[:,'log2FoldChange_x']) != np.sign(merged.loc[:,'log2FoldChange_y'])]
        # theoretically there should be no output here because i == j means we are comparing the same files
        if i == j:
            print(diff)
        # count the number of tfs that appear in both compared lists and change in the same or different directions
        same_gene[i,j] = np.sum(same['_merge'] == 'both')
        diff_gene[i,j] = np.sum(diff['_merge'] == 'both')
        merged = deTFs[i][["Ensembl ID", "gene_name", "log2FoldChange", "padj"]].merge(deTFs[j][["Ensembl ID", "gene_name", "log2FoldChange", "padj"]], on = "Ensembl ID", indicator = True, how = "outer")
        merged.to_csv(os.path.expandvars("$SCRATCH/compareCond/" + outNames[i] + "_" + outNames[j] + "_tfs.csv"))
        same = merged.loc[np.sign(merged.loc[:,'log2FoldChange_x']) == np.sign(merged.loc[:,'log2FoldChange_y'])]
        diff = merged.loc[np.sign(merged.loc[:,'log2FoldChange_x']) != np.sign(merged.loc[:,'log2FoldChange_y'])]
        # count the number of tfs that appear in both compared lists and change in the same or different directions
        same_tf[i,j] = np.sum(same['_merge'] == 'both')
        diff_tf[i,j] = np.sum(diff['_merge'] == 'both')
np.savetxt(os.path.expandvars("$SCRATCH/compareCond/same_gene.csv"), same_gene)
np.savetxt(os.path.expandvars("$SCRATCH/compareCond/diff_gene.csv"), diff_gene)
np.savetxt(os.path.expandvars("$SCRATCH/compareCond/same_tf.csv"), same_tf)
np.savetxt(os.path.expandvars("$SCRATCH/compareCond/diff_tf.csv"), diff_tf)

# merged dataframe would normally go here
ss_ci = pd.DataFrame()

# naive attempt at hypergeometric test to assess significance of the overlap observed between
b = ss_ci.shape[0]
s = ss_tfs.shape[0]
c = ci_tfs.shape[0]
t = tf.shape[0]
pval = hypergeom.sf(b-1, t, s, c)

table = pd.DataFrame({"In SSvP": [b, s-b, s, "NA"], "Not in SSvP": [c-b, t-c-s+b, t-s, "NA"], "Total": [c, t-c, t, pval]}, index = ["In CIvP", "Not in CIvP", "Total", "Hypergeometric p-value"])
table.to_csv(os.path.expandvars("$SCRATCH/compareCond/table.csv"))

plt.scatter(merged.loc[merged["_merge"] == "both"]["l2FC_SSvP"], merged.loc[merged["_merge"] == "both"]["l2FC_CIvP"], c = "steelblue", label = 'DE Genes', alpha = 0.5)
coef = np.polyfit(merged.loc[merged["_merge"] == "both"]["l2FC_SSvP"], merged.loc[merged["_merge"] == "both"]["l2FC_CIvP"], 1)
eq2 = np.poly1d(coef)
r2, p2 = pearsonr(merged.loc[merged["_merge"] == "both"]["l2FC_SSvP"], merged.loc[merged["_merge"] == "both"]["l2FC_CIvP"])
plt.plot(merged.loc[merged["_merge"] == "both"]["l2FC_SSvP"], eq2(merged.loc[merged["_merge"] == "both"]["l2FC_SSvP"]), c = 'r', label = f'r = {round(r2, 3)}')

# plot showing very strong correlation between the log 2 fold changes of genes in all conditions
plt.scatter(ss_ci["l2FC_SSvP"], ss_ci["l2FC_CIvP"], c = 'darkorange', label = 'DE TFs', alpha = 0.5)
plt.xlabel("SSvsP Log 2 Fold Change")
plt.ylabel("CIvsP Log Fold Change")
coef = np.polyfit(ss_ci["l2FC_SSvP"], ss_ci["l2FC_CIvP"], 1)
eq1 = np.poly1d(coef)
r1, p1 = pearsonr(ss_ci["l2FC_SSvP"], ss_ci["l2FC_CIvP"])
plt.plot(ss_ci["l2FC_SSvP"], eq1(ss_ci["l2FC_SSvP"]), c = 'darkorange', label = f'r = {round(r1, 3)}')
plt.legend()
plt.savefig(os.path.expandvars("$SCRATCH/compareCond/corr.png"))