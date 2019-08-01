from rdkit.Chem import AllChem
from rdkit import Chem

ms = []


def ClusterFps(fps, cutoff=0.6):
    from rdkit import DataStructs
    from rdkit.ML.Cluster import Butina
    
    # first generate the distance matrix:
    dists = []
    nfps = len(fps)
    for i in range(1, nfps):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        dists.extend([1 - x for x in sims])
    
    # now cluster the data:
    cs = Butina.ClusterData(dists, nfps, cutoff, isDistData=True)
    return cs


suppl = Chem.SDMolSupplier('/Users/matina/Desktop/tomdata/Master4A.sdf')

for mol in suppl:
    ms.append(mol)
fps = [AllChem.GetMorganFingerprintAsBitVect(x, 2, 2048) for x in ms]  # 2 is radius and 1024b is nBits
clusters = ClusterFps(fps, cutoff=0.5)
file=open('/Users/matina/Desktop/4Aclustering.txt','w+')
print(fps)
print(clusters)
for t in clusters:
    line = ' '.join(str(x) for x in t)
    file.write(line + '\n')
file.close()
print (len(clusters))
