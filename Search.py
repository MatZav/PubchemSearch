from __future__ import print_function
from rdkit import Chem
import os
from progressbar import ProgressBar
pbar=ProgressBar()
matches = []
directory = '/Users/matina/Desktop/'
patt = Chem.MolFromSmarts('NC(N****NC=O)=O')
for file in pbar(os.listdir(directory)):
	filename = os.fsdecode(file)
	if filename.endswith(".sdf"):
		f = os.path.join(directory,filename)
		suppl= Chem.SDMolSupplier(f)
		for mol in suppl:
			if mol is None: continue
			if mol.HasSubstructMatch(patt):
				matches.append(mol)
		w = Chem.SDWriter(Users/matina/Desktop/datasmarts4c.sdf')
		for m in matches: w.write(m)
		print(filename)
		



