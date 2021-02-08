import matplotlib.pyplot as plt
import pandas as pd
import os
from subprocess import call
import numpy as np

atacCsv = os.path.expandvars('$SCRATCH/atac/corr/atacAnno.csv')
atac = os.path.expandvars('$SCRATCH/atac/atacAnno.bed')
pd.read_csv(atacCsv).to_csv(atac, header = None, index = False, sep = '\t')

merged_dir = os.path.expandvars('$SCRATCH/atac/')

# sort is necessary for bedtools merge to work
# merge is to collapse all cCRE annotations for a given DA region into a single line
""" call('sort -k1,1 -k2,2n ' + atac + ' > ' + os.path.join(merged_dir, 'sorted_anno.bed'), shell = True)
call('bedtools merge -i ' + os.path.join(merged_dir, 'sorted_anno.bed') + ' -c 11,13,14 -o mean,mean,distinct > ' + os.path.join(merged_dir, 'clean_anno.bed'), shell = True)

call('sort -k1,1 -k2,2n ' + os.path.join(merged_dir, 'merged_with_encode_ccre.bed') + ' > ' + os.path.join(merged_dir, 'sorted_encode_ccre.bed'), shell = True)
call('bedtools merge -i ' + os.path.join(merged_dir, 'sorted_encode_ccre.bed') + ' -c 9,11,17 -o mean,mean,distinct -delim \",\" > ' + os.path.join(merged_dir, 'clean_encode_ccre.bed'), shell = True)

call('sort -k1,1 -k2,2n ' + os.path.join(merged_dir, 'merged_with_roadmap_chromHMM.bed') + ' > ' + os.path.join(merged_dir, 'sorted_roadmap_chromHMM.bed'), shell = True)
call('bedtools merge -i ' + os.path.join(merged_dir, 'sorted_roadmap_chromHMM.bed') + ' -c 9,11,15 -o mean,mean,distinct -delim \",\" > ' + os.path.join(merged_dir, 'clean_roadmap_chromHMM.bed'), shell = True)

call('sort -k1,1 -k2,2n ' + os.path.join(merged_dir, 'merged_with_roadmap_DHS.bed') + ' > ' + os.path.join(merged_dir, 'sorted_roadmap_DHS.bed'), shell = True)
call('bedtools merge -i ' + os.path.join(merged_dir, 'sorted_roadmap_DHS.bed') + ' -c 9,11,16 -o mean,mean,distinct -delim \",\" > ' + os.path.join(merged_dir, 'clean_roadmap_DHS.bed'), shell = True)
 """

# Drop duplicate annotations (e.g. dELS, dELS) by converting to set
# need to sort list when converting back because sets are UNORDERED
merged_anno = pd.read_csv(os.path.join(merged_dir, 'clean_anno.bed'), sep = '\t', header = None)
unique_anno = [' & '.join(sorted(list(set(i.split(';'))))) for i in merged_anno.iloc[:,5]]
# consolidate all introns and exons into two distinct categories
for count, word in enumerate(unique_anno):
    if word.split(' ')[0] == 'Intron':
        unique_anno[count] = 'Intron'
    elif word.split(' ')[0] == 'Exon':
        unique_anno[count] = 'Exon'
merged_anno.iloc[:,5] = unique_anno
up_anno = merged_anno.loc[merged_anno[3]>=0,]
down_anno = merged_anno.loc[merged_anno[3]<0,]
uanno = pd.Series(up_anno.iloc[:,5])
danno = pd.Series(down_anno.iloc[:,5])
u_freq_anno = uanno.value_counts()/len(uanno)*100
d_freq_anno = danno.value_counts()/len(danno)*100

merged_ccre = pd.read_csv(os.path.join(merged_dir, 'clean_encode_ccre.bed'), sep = '\t', header = None)
unique_ccre = [' & '.join(sorted(list(set(i.split(','))))) for i in merged_ccre.iloc[:,5]]
merged_ccre.iloc[:,5] = unique_ccre
up_ccre = merged_ccre.loc[merged_ccre[3]>=0,]
down_ccre = merged_ccre.loc[merged_ccre[3]<0,]
uccre = pd.Series(up_ccre.iloc[:,5])
dccre = pd.Series(down_ccre.iloc[:,5])
u_freq_ccre = uccre.value_counts()/len(uccre)*100
d_freq_ccre = dccre.value_counts()/len(dccre)*100
u_freq_ccre = u_freq_ccre.rename({'.': 'No Match'})
d_freq_ccre = d_freq_ccre.rename({'.': 'No Match'})

merged_chromHMM = pd.read_csv(os.path.join(merged_dir, 'clean_roadmap_chromHMM.bed'), sep = '\t', header = None)
unique_anno = [', '.join(sorted(list(set(i.split(','))))) for i in merged_chromHMM.iloc[:,5]]
merged_chromHMM.iloc[:,5] = unique_anno
up_chromHMM = merged_chromHMM.loc[merged_chromHMM[3]>=0,]
down_chromHMM = merged_chromHMM.loc[merged_chromHMM[3]<0,]
uchromHMM = pd.Series(up_chromHMM.iloc[:,5])
dchromHMM = pd.Series(down_chromHMM.iloc[:,5])
u_freq_chromHMM = uchromHMM.value_counts()/len(uchromHMM)*100
d_freq_chromHMM = dchromHMM.value_counts()/len(dchromHMM)*100

merged_dhs = pd.read_csv(os.path.join(merged_dir, 'clean_roadmap_DHS.bed'), sep = '\t', header = None)
unique_anno = [', '.join(sorted(list(set(i.split(','))))) for i in merged_dhs.iloc[:,5]]
merged_dhs.iloc[:,5] = unique_anno
up_dhs = merged_dhs.loc[merged_dhs[3]>=0,]
down_dhs = merged_dhs.loc[merged_dhs[3]<0,]
udhs = pd.Series(up_dhs.iloc[:,5])
ddhs = pd.Series(down_dhs.iloc[:,5])
u_freq_dhs = udhs.value_counts()/len(udhs)*100
d_freq_dhs = ddhs.value_counts()/len(ddhs)*100

fig, ax1 = plt.subplots()
freq_anno = pd.concat([u_freq_anno, d_freq_anno], axis = 1)
for i in range(freq_anno.shape[0]):
    if i + 1 < 7:
        ax1.barh((1,0), freq_anno.iloc[i,:].tolist(), left = freq_anno[i+1:].sum(axis = 0).tolist(), label = freq_anno.index[i])
    else:
        ax1.barh((1,0), freq_anno[i:].sum(axis = 0).tolist(), label = 'Other')
        break
ax1.set_xlabel('%')
lgd = ax1.legend(loc='center left', bbox_to_anchor=(1,0.5))
ax1.set_yticks([1,0])
ax1.set_yticklabels(['Up', 'Down'])
fig.set_size_inches(5, 3)
fig.savefig(merged_dir + "figures/anno.png", bbox_extra_artists=[lgd], bbox_inches='tight')

fig, ax1 = plt.subplots()
freq_ccre = pd.concat([u_freq_ccre, d_freq_ccre], axis = 1)
for i in range(freq_ccre.shape[0]):
    if i + 1 < 7:
        ax1.barh((1,0), freq_ccre.iloc[i,:].tolist(), left = freq_ccre[i+1:].sum(axis = 0).tolist(), label = freq_ccre.index[i])
    else:
        ax1.barh((1,0), freq_ccre[i:].sum(axis = 0).tolist(), label = 'Other')
        break
ax1.set_xlabel('%')
lgd = ax1.legend(loc='center left', bbox_to_anchor=(1,0.5))
ax1.set_yticks([1,0])
ax1.set_yticklabels(['Up', 'Down'])
fig.set_size_inches(5, 3)
fig.savefig(merged_dir + "figures/ccre.png", bbox_extra_artists=[lgd], bbox_inches='tight')


def match(old):
    old = old.split(', ')
    for i, state in enumerate(old):
        old[i] = ''.join(statesKey.loc[statesKey["MNEMONIC"]==state, "DESCRIPTION"].tolist())
    return ' & '.join(old)

# covert esoteric chromatin state names to more parsable format
statesKey = pd.read_csv(merged_dir + 'states.txt', sep = '\t')
statesKey["MNEMONIC"] = [str(statesKey.iloc[i,0]) + '_' + statesKey.iloc[i,1] for i in statesKey.index]
u_freq_chromHMM = u_freq_chromHMM.rename(match)
d_freq_chromHMM = d_freq_chromHMM.rename(match)
fig, ax1 = plt.subplots()
freq_chromHMM = pd.concat([u_freq_chromHMM, d_freq_chromHMM], axis = 1)
for i in range(freq_chromHMM.shape[0]):
    if i + 1 < 7:
        ax1.barh((1,0), freq_chromHMM.iloc[i,:].tolist(), left = freq_chromHMM[i+1:].sum(axis = 0).tolist(), label = freq_chromHMM.index[i])
    else:
        ax1.barh((1,0), freq_chromHMM[i:].sum(axis = 0).tolist(), label = 'Other')
        break
ax1.set_xlabel('%')
lgd = ax1.legend(loc='center left', bbox_to_anchor=(1,0.5))
ax1.set_yticks([1,0])
ax1.set_yticklabels(['Up', 'Down'])
fig.set_size_inches(5, 3)
fig.savefig(merged_dir + "figures/chromHMM.png", bbox_extra_artists=[lgd], bbox_inches='tight')

print(len(uanno))
print(len(danno))