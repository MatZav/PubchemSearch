from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Lipinski
from rdkit.Chem import Descriptors
w4a = Chem.SDWriter('/Users/matina/Desktop/MasterReduce.sdf')
suppl = Chem.SDMolSupplier('/Users/matina/Desktop/Master.sdf')
len(suppl)

for mol in suppl:
    if mol is None: continue
    a = Lipinski.NumHDonors(mol)
    b = Descriptors.MolWt(mol)
    if a>2 and b<1000:
        w4a.write(mol)