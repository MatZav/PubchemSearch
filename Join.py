import csv
import os
from rdkit import Chem
from rdkit.Chem import AllChem

C99 = {'H': '[H]%99',
       'Me': 'C%99',
       '3,5-(CF3)2C6H3': 'C1%99=CC(C(F)(F)F)=CC(C(F)(F)F)=C1',
       'Br': '[Br]%99'}
C98 = {'H': 'H%98',
       'Me': 'C%98',
       '3,5-(CF3)2C6H3': 'C1%98=CC(C(F)(F)F)=CC(C(F)(F)F)=C1',
       'Br': 'Br%98'}
C97 = {'Ph': 'C1%97=CC=CC=C1',
       'Cy': 'C1%97CCCCC1',
       '3,5-(CF3)2C6H3': 'C1%97=CC(C(F)(F)F)=CC(C(F)(F)F)=C1',
       ' 4-(NO2)C6H5': 'C1%97=CC=C([N+]([O-])=O)C=C1',
       ' 1-naphthyl': 'C1%97=C(C=CC=C2)C2=CC=C1',
       '2-(CF3)C6H4': 'C1%97=C(C(F)(F)F)C=CC=C1',
       '2,6-(Cl)2C6H3': 'C1%97=C(Cl)C=CC=C1Cl',
       ' Ts': 'S%97(=O)(C1=CC=C(C)C=C1)=O',
       'C6F5': 'C1%97=C(F)C(F)=C(F)C(F)=C1F',
       '3,5-[3,5-(CF3)2C6H3]2C6H3': 'C1%97=CC(C2=CC(C(F)(F)F)=CC(C(F)(F)F)=C2)=CC(C3=CC(C(F)(F)F)=CC(C(F)(F)F)=C3)=C1',
       '3,5-[3,5-(C6H5)2C6H3]2C6H3': 'C1%97=CC(C2=CC(C3=CC=CC=C3)=CC(C4=CC=CC=C4)=C2)=CC(C5=CC(C6=CC=CC=C6)=CC(C7=CC=CC=C7)=C5)=C1',
       '3,5-[3,5-(CH3)2C6H3]2C6H3': 'C1%97=CC(C2=CC(C)=CC(C)=C2)=CC(C3=CC(C)=CC(C)=C3)=C1',
       '3,5-[2,5-(CF3)2C6H3]2C6H3': 'C1%97=CC(C2=CC(C(F)(F)F)=CC=C2C(F)(F)F)=CC(C3=C(C(F)(F)F)C=CC(C(F)(F)F)=C3)=C1',
       '4-(SF5)C6H5': 'C1%97=CC=C(S(F)(F)(F)(F)F)C=C1',
       '3,5-(CH3)2C6H5': 'C1%97=CC(C)=CC(C)=C1',
       '4-CH3-3,5-[3,5-(CF3)2C6H3]2C6H3': 'C1%97=CC(C2=CC(C(F)(F)F)=CC(C(F)(F)F)=C2)=C(C)C(C3=CC(C(F)(F)F)=CC(C(F)(F)F)=C3)=C1',
       ' 4-F-3,5-[3,5-(CF3)2C6H3]2C6H3': 'C1%97=CC(C2=CC(C(F)(F)F)=CC(C(F)(F)F)=C2)=C(F)C(C3=CC(C(F)(F)F)=CC(C(F)(F)F)=C3)=C1',
       '3,5-[4-Me-3,5-(CF3)2C6H3]2C6H3': 'C1%97=CC(C2=CC(C(F)(F)F)=C(C)C(C(F)(F)F)=C2)=CC(C3=CC(C(F)(F)F)=C(C)C(C(F)(F)F)=C3)=C1',
       '3,5-[4-F-3,5-(CF3)2C6H3]2C6H3': 'C1%97=CC(C2=CC(C(F)(F)F)=C(F)C(C(F)(F)F)=C2)=CC(C3=CC(C(F)(F)F)=C(F)C(C(F)(F)F)=C3)=C1',
       '3,5-[3,5-(i-C3F7)2C6H3]2C6H3': 'C1%97=CC(C2=CC(C(F)(C(F)(F)F)C(F)(F)F)=CC(C(C(F)(F)F)(C(F)(F)F)F)=C2)=CC(C3=CC(C(F)(C(F)(F)F)C(F)(F)F)=CC(C(F)(C(F)(F)F)C(F)(F)F)=C3)=C1',
       '3,5-[3,5-(SF5)2C6H3]2C6H3': 'C1%97=CC(C2=CC([SF5])=CC(C)=C2)=CC(C3=CC([SF5])=CC([SF5])=C3)=C1',
       '3,5-(SF5)2C6H': 'C1%97=CC(S(F)(F)(F)(F)F)=CC(S(F)(F)(F)(F)F)=C1',
       '3,5-(i-C3F7)2C6H3': 'C1%97=CC(C(F)(C(F)(F)F)C(F)(F)F)=CC(C(F)(C(F)(F)F)C(F)(F)F)=C1'}

C96 = {'Ph': 'C1%96=CC=CC=C1',
       'Cy': 'C1%96CCCCC1',
       '3,5-(CF3)2C6H3': 'C1%96=CC(C(F)(F)F)=CC(C(F)(F)F)=C1',
       ' 4-(NO2)C6H5': 'C1%96=CC=C([N+]([O-])=O)C=C1',
       ' 1-naphthyl': 'C1%96=C(C=CC=C2)C2=CC=C1',
       '2-(CF3)C6H4': 'C1%96=C(C(F)(F)F)C=CC=C1',
       '2,6-(Cl)2C6H3': 'C1%96=C(Cl)C=CC=C1Cl',
       ' Ts': 'S%96(=O)(C1=CC=C(C)C=C1)=O',
       'C6F5': 'C1%96=C(F)C(F)=C(F)C(F)=C1F',
       '3,5-[3,5-(CF3)2C6H3]2C6H3': 'C1%96=CC(C2=CC(C(F)(F)F)=CC(C(F)(F)F)=C2)=CC(C3=CC(C(F)(F)F)=CC(C(F)(F)F)=C3)=C1',
       '3,5-[3,5-(C6H5)2C6H3]2C6H3': 'C1%96=CC(C2=CC(C3=CC=CC=C3)=CC(C4=CC=CC=C4)=C2)=CC(C5=CC(C6=CC=CC=C6)=CC(C7=CC=CC=C7)=C5)=C1',
       '3,5-[3,5-(CH3)2C6H3]2C6H3': 'C1%96=CC(C2=CC(C)=CC(C)=C2)=CC(C3=CC(C)=CC(C)=C3)=C1',
       '3,5-[2,5-(CF3)2C6H3]2C6H3': 'C1%96=CC(C2=CC(C(F)(F)F)=CC=C2C(F)(F)F)=CC(C3=C(C(F)(F)F)C=CC(C(F)(F)F)=C3)=C1',
       '4-(SF5)C6H5': 'C1%96=CC=C(S(F)(F)(F)(F)F)C=C1',
       '3,5-(CH3)2C6H5': 'C1%96=CC(C)=CC(C)=C1',
       '4-CH3-3,5-[3,5-(CF3)2C6H3]2C6H3': 'C1%96=CC(C2=CC(C(F)(F)F)=CC(C(F)(F)F)=C2)=C(C)C(C3=CC(C(F)(F)F)=CC(C(F)(F)F)=C3)=C1',
       ' 4-F-3,5-[3,5-(CF3)2C6H3]2C6H3': 'C1%96=CC(C2=CC(C(F)(F)F)=CC(C(F)(F)F)=C2)=C(F)C(C3=CC(C(F)(F)F)=CC(C(F)(F)F)=C3)=C1',
       '3,5-[4-Me-3,5-(CF3)2C6H3]2C6H3': 'C1%96=CC(C2=CC(C(F)(F)F)=C(C)C(C(F)(F)F)=C2)=CC(C3=CC(C(F)(F)F)=C(C)C(C(F)(F)F)=C3)=C1',
       '3,5-[4-F-3,5-(CF3)2C6H3]2C6H3': 'C1%96=CC(C2=CC(C(F)(F)F)=C(F)C(C(F)(F)F)=C2)=CC(C3=CC(C(F)(F)F)=C(F)C(C(F)(F)F)=C3)=C1',
       '3,5-[3,5-(i-C3F7)2C6H3]2C6H3': 'C1%96=CC(C2=CC(C(F)(C(F)(F)F)C(F)(F)F)=CC(C(C(F)(F)F)(C(F)(F)F)F)=C2)=CC(C3=CC(C(F)(C(F)(F)F)C(F)(F)F)=CC(C(F)(C(F)(F)F)C(F)(F)F)=C3)=C1',
       '3,5-[3,5-(SF5)2C6H3]2C6H3': 'C1%96=CC(C2=CC([SF5])=CC(C)=C2)=CC(C3=CC([SF5])=CC([SF5])=C3)=C1',
       '3,5-(SF5)2C6H': 'C1%96=CC(S(F)(F)(F)(F)F)=CC(S(F)(F)(F)(F)F)=C1',
       '3,5-(i-C3F7)2C6H3': 'C1%96=CC(C(F)(C(F)(F)F)C(F)(F)F)=CC(C(F)(C(F)(F)F)C(F)(F)F)=C1'}

C95 = {'H': 'H%95',
       'Br': 'Br%95',
       'I': 'I%95',
       ' 2-naphthyl': 'C1%95=CC(C=CC=C2)=C2C=C1',
       ' Ph': 'C1%95=CC=CC=C1',
       ' 2-(CH3)C6H4': 'C1%95=C(C)C=CC=C1',
       ' 4-(CH3)C6H4': 'C1%95=CC(C(F)(F)F)=CC(C(F)(F)F)=C1',
       '3,5-(CF3)2C6H3': 'C1%95=CC(C(F)(F)F)=CC(C(F)(F)F)=C1',
       '4-(CF3)C6H4': 'C1%95=CC=C(C(F)(F)F)C=C1',
       ' 4-(C6H5)C6H4': 'C1%95=CC=C(C2=CC=CC=C2)C=C1',
       '4-(tBu)C6H4': 'C1%95=CC=C(C(C)(C)C)C=C1',
       ' 2,5-(CH3)2C6H3': 'C1%95=C(C)C=CC(C)=C1',
       '(4-(1-naphthyl)C6H4': 'C1%95=CC=C(C2=C(C=CC=C3)C3=CC=C2)C=C1',
       '(1-(2-naphthyl)C6H4': 'C1%95=C(C2=CC=CC3=C2C=CC=C3)C=CC=C1',
       '2-(OH)C6H4': 'C1%95=C(O)C=CC=C1',
       '4-(4-(Et)C6H4)C6H4': 'C1%95=CC=C(C2=CC=C(CC)C=C2)C=C1',
       ' 4-(OMe)C6H5': 'C1%95=CC=C(OC)C=C1',
       '3,5-(CH3)2C6H3': 'C1%95=CC=C(OC)C=C1'}
C94 = {'H': 'H%94',
       'Br': 'Br%94',
       'I': 'I%94',
       ' 2-naphthyl': 'C1%94=CC(C=CC=C2)=C2C=C1',
       ' Ph': 'C1%94=CC=CC=C1',
       ' 2-(CH3)C6H4': 'C1%94=C(C)C=CC=C1',
       ' 4-(CH3)C6H4': 'C1%94=CC(C(F)(F)F)=CC(C(F)(F)F)=C1',
       '3,5-(CF3)2C6H3': 'C1%94=CC(C(F)(F)F)=CC(C(F)(F)F)=C1',
       '4-(CF3)C6H4': 'C1%94=CC=C(C(F)(F)F)C=C1',
       ' 4-(C6H5)C6H4': 'C1%94=CC=C(C2=CC=CC=C2)C=C1',
       '4-(tBu)C6H4': 'C1%94=CC=C(C(C)(C)C)C=C1',
       ' 2,5-(CH3)2C6H3': 'C1%94=C(C)C=CC(C)=C1',
       '(4-(1-naphthyl)C6H4': 'C1%94=CC=C(C2=C(C=CC=C3)C3=CC=C2)C=C1',
       '(1-(2-naphthyl)C6H4': 'C1%94=C(C2=CC=CC3=C2C=CC=C3)C=CC=C1',
       '2-(OH)C6H4': 'C1%94=C(O)C=CC=C1',
       '4-(4-(Et)C6H4)C6H4': 'C1%94=CC=C(C2=CC=C(CC)C=C2)C=C1',
       ' 4-(OMe)C6H5': 'C1%94=CC=C(OC)C=C1',
       '3,5-(CH3)2C6H3': 'C1%94=CC=C(OC)C=C1'}
C93 = {'H': 'H%93', 'Me': 'C%93', 'Et': 'CC%93', 'Pr': 'CCC%93', 'iPr': 'C(C)C%93'}

core1 = 'O=C(N%97)N%93C1=[C@]([C@@]2=C(C=CC%99=C3)C3=CC%95=C2NC(N%96)=O)C4=CC=C%98C=C4C=C1%94'
core2 = 'O=C(N%97)N%93C1=[C@]([C@@]2=C(CCC%99C3)C3=CC%95=C2NC(N%96)=O)C(CCC%98C4)=C4C=C1%94'


def make_linker_smilesv2(c99, c98, c97, c96, c95, c94,c93,core):
    A=C99[c99]
    B=C98[c98]
    C=C97[c97]
    D=C96[c96]
    E=C95[c95]
    F=C94[c94]
    G=C93[c93]
    H=core
    return ( A +'.' + B+'.'+C+'.'+D+'.'+E+'.'+F+'.'+G +'.' + H)


l=[]

for file in os.listdir('Users/matina/Desktop'):
    filename = os.fsdecode(file)
    if filename.endswith(".csv"):
        l.append(filename)

b = make_linker_smilesv2('Br','Me','Ph','Ph',' Ph',' Ph','Me',core2)
suppl = Chem.SDMolSupplier('/Users/matina/Desktop/A.smi')
mol = [x for x in suppl]
mol = Chem.MolFromSmiles(b)
print(Chem.MolToSmiles(mol))