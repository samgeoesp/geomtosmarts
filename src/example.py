from geomtosmarts import geo

smiles, smarts = geo('../examples/diels_alder.out', r=True, k=True, n=True)

print(f'Reaction SMARTS:\n{smarts}')
print('---')
print(f'Reaction SMILES:\n{smiles}')


"""
Command Line Usage:

python geomtosmarts.py -f ../examples/diels_alder.out -r 

"""