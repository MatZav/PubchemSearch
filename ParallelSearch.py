from __future__ import print_function
from rdkit import Chem
import os
from progressbar import ProgressBar
from joblib import Parallel, delayed
import multiprocessing
from multiprocessing import  Process, Queue

pbar=ProgressBar()
matches = []
directory = 'E:\Data'
patt = Chem.MolFromSmarts('[N;H1]C=O.[N;H1]C([N;H1])=O')
w = 'E:\Search3'
l=[]

for file in pbar(os.listdir(directory)):
    filename = os.fsdecode(file)
    if filename.endswith(".sdf"):
        l.append(filename)

num_cores = multiprocessing.cpu_count()
print(num_cores)
lock = multiprocessing.Lock()

def Search(i,w,directory):
    b = os.path.join(w,i)
    j=os.path.join(directory,i)
    suppl= Chem.SDMolSupplier(j)
    g = Chem.SDWriter(b)
    for mol in suppl:
        if mol is None: continue
        if mol.HasSubstructMatch(patt):
            g.write(mol)
    return

results = Parallel(n_jobs=20)(delayed(Search)(i,w,directory) for i in l)
k=[]
filename=[]
w = 'E:\Search3'

for file in os.listdir('E:\Search3'):
    filename = os.fsdecode(file)
    if filename.endswith(".sdf"):
        k.append(filename)

for file_name in k:
    file_name = os.path.join(w,file_name)
    # This needs to be done *inside the loop*
    f= open(file_name, 'r')
    lst = []
    for line in f:
       lst.append(line)
    f.close()

    f=open('E:\\Search3\\Master.sdf', 'a')

    for line in lst:
        f.write(line)
    f.close()

suppl= Chem.SDMolSupplier('E:\\Search3\\Master.sdf')