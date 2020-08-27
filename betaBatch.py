import subprocess
import shlex
import os
import re
from multiprocessing import Pool

peakdir = os.path.expandvars('$SCRATCH/tfchip/goodBeds/')
outdir = os.path.expandvars('$SCRATCH/output/betaMinus/')

def beta(name, root):
    args = shlex.split('BETA minus -p ' + name + ' -g hg38 -n ' + os.path.dirname(os.path.join(root, name)) + '_' + re.sub('[^0-9]', '', name) + ' -o ' + outdir)
    print(subprocess.call(args, universal_newlines=True, cwd=root))

def unpack(args):
    beta(*args)

names = []
roots = []
for root, dirs, files in os.walk(peakdir):
    for name in files:
        names.append(name)
        roots.append(root)
names.remove("nodata.csv")
names.remove("key.csv")
roots.remove(peakdir)
roots.remove(peakdir)
p = Pool()
p.map(unpack, zip(names, roots))