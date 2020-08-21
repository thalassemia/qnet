import os
import glob
import pandas
from shutil import copy

os.chdir('/mnt/c/Users/Sean/Desktop/Research/')
def iterable(obj):
    try:
        iter(obj)
    except Exception:
        return False
    else:
        return True

fail = []

sigTFs = pandas.read_csv('SSvsP_SigTFs.csv', index_col='gene_name')
index = pandas.read_csv('Human TF ChIP-Seq Index.txt', sep='\t', index_col='Factor')
os.chdir('/mnt/c/Users/Sean/Desktop/Research/Human TF ChIP-Seq/')
for geneName in sigTFs.index:
    if (index.index == geneName).any():
        nums = index.at[geneName, 'DCid']
        if(iterable(nums)):
            for num in nums:
                for file in glob.glob(f'{num}*.bed'):
                    os.makedirs(os.path.dirname(f'/home/sean/peaks/{geneName}/'), exist_ok=True)
                    copy(file, f'/home/sean/peaks/{geneName}/')
        else:
            for file in glob.glob(f'{nums}*.bed'):
                os.makedirs(os.path.dirname(f'/home/sean/peaks/{geneName}/'), exist_ok=True)
                copy(file, f'/home/sean/peaks/{geneName}/')
    else:
        fail.append(geneName)
NoData = pandas.DataFrame(fail)
pandas.to_csv(NoData.T)