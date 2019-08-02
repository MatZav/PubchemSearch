import os
import sys
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.spatial import distance
import numpy as np
import progressbar
from joblib import Parallel, delayed
import multiprocessing
from multiprocessing import  Process, Queue

x,y,z=[],[],[]
flag,count=0,0
numConfs=10
coordinatesMatrix=[]
distancesMatrix=[]
value3A=[]
value4A=[]
matches3A=[]
matches4A=[]
l=[]

directory = 'E:\Distances'
A_3= 'E:\\Distances\\3A'
A_4= 'E:\\Distances\\4A'

for file in os.listdir(directory):
    filename = os.fsdecode(file)
    if filename.endswith(".sdf"):
        l.append(filename)
i=l[1]
w=os.path.join(directory,i)
A3 = os.path.join(A_3,i)
A4 = os.path.join(A_4,i)

num_cores = multiprocessing.cpu_count()
print(num_cores)
lock = multiprocessing.Lock()


def Distance(i, directory, A_3, A_4, count):
    coordinatesMatrix = []
    distancesMatrix = []
    value3A = []
    value4A = []
    matches3A = []
    matches4A = []
    w = os.path.join(directory, i)
    A3 = os.path.join(A_3, i)
    A4 = os.path.join(A_4, i)
    suppl = Chem.SDMolSupplier(w)
    w3A = Chem.SDWriter(A3)
    w4A = Chem.SDWriter(A4)
    for mol in suppl:
        patt = Chem.MolFromSmarts('NC(N)=O')
        num = mol.GetSubstructMatches(patt)
        h = len(num)
        m3 = Chem.AddHs(mol)
        cids = AllChem.EmbedMultipleConfs(m3, numConfs)
        for i in cids:
            a = Chem.MolToMolBlock(m3, confId=i)
            for item in a.split('\n'):
                if 'N' in item:
                    if 'END' in item:
                        continue
                    else:
                        cord = item.strip(' ')

                        seqx = (cord[0], cord[1], cord[2], cord[3], cord[4], cord[5], cord[6])
                        jx = ''.join(seqx)
                        x.append(jx)

                        seqy = (cord[9], cord[10], cord[11], cord[12], cord[13], cord[14], cord[15])
                        jy = ''.join(seqy)
                        y.append(jy)

                        seqz = (cord[19], cord[20], cord[21], cord[22], cord[23], cord[24], cord[25])
                        jz = ''.join(seqz)
                        z.append(jz)

                        count = count + 1
        # count=0
        # sys.stdout = open("out.txt", "w")
        # sys.stdout.close()
        cc = int(count / numConfs)

        float_x = [float(ix) for ix in x]
        float_y = [float(iy) for iy in y]
        float_z = [float(iz) for iz in z]
        myArray = np.array((float_x, float_y, float_z))
        myArray.astype(float)

        for i in range(count):
            a = myArray[0][i]
            b = myArray[1][i]
            c = myArray[2][i]
            coordinates = (a, b, c)
            coordinatesMatrix.append(coordinates)

        i = 0
        w = cc - 1
        w2 = cc - 2
        w3 = cc - 1

        for p in range(numConfs):

            while w > 0:
                dst = distance.euclidean(coordinatesMatrix[i], coordinatesMatrix[i + w])
                distancesMatrix.append(dst)
                w = w - 1

            while w2 > 0 and w3 > 1:
                dst = distance.euclidean(coordinatesMatrix[i + w2], coordinatesMatrix[i + w3])
                distancesMatrix.append(dst)
                w2 = w2 - 1
                if w2 == 0:
                    break
                else:
                    dst = distance.euclidean(coordinatesMatrix[i + w2], coordinatesMatrix[i + w3])
                    distancesMatrix.append(dst)
                    w3 = w3 - 1

            w = cc - 1
            w2 = cc - 2
            w3 = cc - 1
            i = i + cc

        for item in distancesMatrix:
            if item < 3:
                value3A.append(item)
        if cc == 0:
            continue
        else:
            j = (count / cc) * h
            if len(value3A) > j:
                w3A.write(mol)

        for item in distancesMatrix:
            if item < 4:
                value4A.append(item)
        if cc == 0:
            continue
        else:
            j = (count / cc) * h
            if len(value4A) > j:
                w4A.write(mol)

        del distancesMatrix[:]
        # del coordinatesMatrix
        del coordinatesMatrix[:]
        # del distancesMatrix
        del x[:]
        # del x
        del y[:]
        # del y
        del z[:]
        # del z
        count = 0
        a = []
        j = []
        value3A = []
        value4A = []


    return

results = Parallel(n_jobs=20)(delayed(Distance)(i,directory,A_3,A_4,count) for i in l)

k=[]
w='E:\\Distances\\4A'

for file in os.listdir(w):
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

    f=open('/Users/matina/Desktop/Master4A.sdf', 'a')

    for line in lst:
        f.write(line)
    f.close()

suppl = Chem.SDMolSupplier('/Users/matina/Desktop/Master4A.sdf')
