import subprocess
import shlex
import os
import re
import concurrent.futures

peakdir = '/home/sean/peaks/'

def beta(name, root):

    args = shlex.split('BETA minus -p ' + name + ' -k BSF --gname2 -g hg38 --da 500 -n ' + os.path.dirname(os.path.join(root, name)) + '_' + re.sub('[^0-9]', '', name) + ' -o /home/sean/targets/')
    print(name + subprocess.call(args, universal_newlines=True, cwd=root))


names = []
roots = []
for root, dirs, files in os.walk(peakdir):
    for name in files:
        names.append(name)
        roots.append(root)
with concurrent.futures.ProcessPoolExecutor() as executor:
    executor.map(beta, names, roots)