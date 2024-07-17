"""
Example code for how to use geomtosmarts.py within a python script.
"""

# Import packages.
import os
from geomtosmarts import geo


# Only loading in one structure.
smiles, smarts = geo('../examples/diels_alder.out', r=False, k=False, n=False)
print(f'Reaction SMARTS:\n{smarts}')
print('---')
print(f'Reaction SMILES:\n{smiles}')


# Loading multiple files.
data_dict = {}
for file in os.listdir('../examples/'):
    path = f'../examples/{file}'
    smiles, smarts = geo(path, r=False, k=False, n=False)
    data_dict[file] = {'rxn_smiles':smiles, 'rxn_smarts':smarts}
    
print('---')
print(f'Files passed = {len(data_dict)}')
print('---')
print(f'Data:\n{data_dict}')

"""
Command Line Usage:

python geomtosmarts.py -f ../examples/diels_alder.out -r 

"""