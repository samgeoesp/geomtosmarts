#####################################################################################
#                                                                                   #
#                                                                                   #
#                               geomtosmarts.py                                     #
#                                                                                   #
#        Code to take an optimised Gaussian TS structure and to generate a          #
#        reaction SMARTS.                                                           #
#                                                                                   #
#              Usage:                                                               #
#                     python geomtosmarts.py                                        #
#                                                                                   #
#              Arguments:                                                           #
#                     -r --> Tag to save .png of reaction SMARTS.                   #
#                     -k --> Tag to keep the intermediate .xyz files.               #
#                     -n --> Tag to keep the file name with saved outputs.          #
#                                                                                   #
#                                                                                   #
#####################################################################################
#                                                                                   #
#                         Written by Samuel G Espley                                #
#                                                                                   #
#####################################################################################


# Imports
import os
import argparse
import numpy as np
import xyz_py as xyzp
from rdkit import Chem
from typing import Tuple
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem import rdDetermineBonds
from molml.utils import get_connections
from molml.constants import BOND_LENGTHS
from cclib.io import ccread # type: ignore

# Add Br length to molml constants.
BOND_LENGTHS['Br'] = {"1":1.3}

def _get_data(file) -> dict:
    """Function to take all data in current working directory and parse
    the appropriate data using cclib.

    Returns:

    data_dict (dictionary): Dictionary with filename (key) and parsed data (value) pairs.
    """
    
    atomic_masses = {1:'H', 6:'C', 7:'N', 8:'O', 9:'F', 14:'Si', 17:'Cl', 35:'Br'}
    data_dict = {}
    d = {}
    # Loop through files and extract information with cclib.
    d[file] = ccread(file)
    for vib in d[file].vibfreqs:
        if vib < 0:
            position = d[file].vibfreqs.tolist().index(vib)
    info_d = {}
    # Pull out coordinates.
    info_d['coords1'] = d[file].atomcoords[-1].astype(float)
    info_d['coords2'] = d[file].atomcoords[-1].astype(float) - d[file].vibdisps[position]
    info_d['coords3'] = d[file].atomcoords[-1].astype(float) + d[file].vibdisps[position]
    # Pull out maximum allowed cooordinates to allow retry.
    info_d['max_coords2'] = d[file].atomcoords[-1].astype(float) - d[file].vibdisps[position]*1.1
    info_d['max_coords3'] = d[file].atomcoords[-1].astype(float) + d[file].vibdisps[position]*1.1
    info_d['disp_vect'] = d[file].vibdisps[position]
    info_d['atomnos'] = d[file].atomnos
    info_d['charge'] = d[file].charge
    info_d['atomcharges'] = d[file].atomcharges['apt']
    atom_symbs = []
    for atom in info_d['atomnos']:
        atom_symbs.append(atomic_masses[atom])
    info_d['atomsymb'] = np.array(atom_symbs)
        
    data_dict[file] = info_d
    return data_dict

def _get_adjacency(data_dict: dict) -> Tuple[dict, dict]:
    """Function to get the adjacency matrix for both displacement coordinates and determine the
    bond forming atoms.

    Arguments:
    data_dict (Dictionary): A dictionary with filename (key) and parsed data (value) pairs.

    Returns:
    data_dict (Dictionary): A dictionary with filename (key) and parsed data (value) pairs.
    react_dict (Dictionary): A dictionary of the structures and the reaction centres.
    """

    # Build empty dictionary.
    react_dict = {}
    # Loop through all structures.
    for structure in data_dict.keys():
        # Get elements and coordinates.
        e_l = list(data_dict[structure]['atomsymb'])
        c1_l = [list(i) for i in data_dict[structure]['coords2']]
        c2_l = [list(i) for i in data_dict[structure]['coords3']]
        # Get the adjacency matrices.
        m1 = xyzp.get_adjacency(e_l, c1_l, 0) # type: ignore
        m2 = xyzp.get_adjacency(e_l, c2_l, 0) # type: ignore
        if np.sum(m1) > np.sum(m2):
            p_c = 'coords2'
            r_c = 'coords3'
        else:
            p_c = 'coords3'
            r_c = 'coords2'
        # Track the changes from one matrix to another.
        changes = []
        for i in range(0, len(m1)):
            if list(m1[i]) != list(m2[i]):
                for index, (fi, se) in enumerate(zip(list(m1[i]), list(m2[i]))):
                    if fi != se:
                        changes.append(f'{i}_{index}')
        # Pull out the first instance and use that. Potentially may have to alter this depending upon the
        # if the structure has unusual behaviour but, this should conistently work.
        if changes == []:
            print(f'No TS found for structure {structure}')
        else:
            react_dict[structure] = changes[0]
            data_dict[structure]['product_coords'] = p_c
            data_dict[structure]['reactant_coords'] = r_c

    return data_dict, react_dict

def _get_connected_atoms(data_dict: dict) -> Tuple[dict, dict]:
        """Function to get the full connectivity of the whole molecule from the perspective of the reacting atom.

        Arguments:
        data_dict (Dictionary): A dictionary with filename (key) and parsed data (value) pairs.

        Returns:
        data_dict (Dictionary): A dictionary with filename (key) and parsed data (value) pairs.
        connect_dict (Dictionary): A dictionary containing the list of atoms that are connected to the reacting atom.
        """
        # Create the distance/vector dictionaries and the reacting atoms for every structure.
        data_dict, react_dict = _get_adjacency(data_dict)
        connect_dict = {}

        # Loop through every structure
        for structure, atoms in react_dict.items():
            connected_atoms = []
            already_checked = []
            # Get the first reacting atom for the current structure and calculate its connectivity storing it in connected_atoms/already_checked
            atom_one = atoms.split('_')[0] 
            connectivity = get_connections(data_dict[structure]['atomsymb'].tolist(), data_dict[structure]['coords1'].tolist())
            for connected_atom in connectivity[int(atom_one)].keys():
                connected_atoms.append(connected_atom)
            connected_atoms.append(int(atom_one))
            already_checked.append(int(atom_one))
            # Now loop through the connected atoms to get the next set of connected atoms checking if its in already_checked to avoid extra computation
            for new_atom in connected_atoms:
                if new_atom not in already_checked:
                    for new_a in connectivity[int(new_atom)].keys():
                        connected_atoms.append(new_a)
                    already_checked.append(new_atom)
            connect_dict[structure] = [*set(connected_atoms)] # Removes duplicate atoms in the list

        return data_dict, connect_dict

def _xyz_creator(data_dict: dict, connect_dict: dict) -> dict:
    """Function to create temporary xyz files for converting to SMILES/SMARTS.

    Arguments:
    data_dict (Dictionary): Dictionary with filename (key) and parsed data (value) pairs.
    connect_dict (Dictionary): A dictionary containing the list of atoms that are connected to the reacting atom.

    Returns:
    rxn_dict (Dictionary): A dictionary containing the information on the reactants/products for a given reaction.
    """
    # Loop through each structure in the connect_dict
    rxn_dict = {}
    for structure in connect_dict.keys():
        # Clean filename.
        filename = os.path.basename(structure)
        name1 = filename.split('.')[0]+'_reactant_1.xyz'
        file_1 = open(name1, 'w')
        file_1.write(f'{len(connect_dict[structure])}\n\n')
        # Loop through all atoms in the connect_dict for the current structure and pull the coordinates from data_dict
        rct1_charge_list = []
        for atom in connect_dict[structure]:
            coords = list(data_dict[structure][data_dict[structure]['reactant_coords']][atom])
            coords = ' '.join([str(float(i)) for i in list(coords)])
            rct1_charge_list.append(list(data_dict[structure]['atomcharges'])[atom])
            coords = f"{str(data_dict[structure]['atomsymb'][atom])} {coords}\n"
            file_1.write(coords)
        file_1.close()
        
        # Here we want to get all the atoms in the TS in a list to remove reactant_1 atoms to leave us with reactant_2 atoms
        name2 = filename.split('.')[0]+'_reactant_2.xyz'
        all_atoms = list(range(0, len(data_dict[structure]['coords1'])))
        second_connect = list(set(all_atoms)- set(connect_dict[structure]))

        # Repeat above but this time for reactant_2
        file_2 = open(name2, 'w')
        file_2.write(f'{len(second_connect)}\n\n')
        rct2_charge_list = []
        for atom in second_connect:
            coords = list(data_dict[structure][data_dict[structure]['reactant_coords']][atom])
            coords = ' '.join([str(float(i)) for i in list(coords)])
            rct2_charge_list.append(list(data_dict[structure]['atomcharges'])[atom])
            coords = f"{str(data_dict[structure]['atomsymb'][atom])} {coords}\n"
            file_2.write(coords)
        file_2.close()

        # Repeat the same for the product
        p_name =  filename.split('.')[0]+'_product.xyz'
        file_3 = open(p_name, 'w')
        file_3.write(f"{len(data_dict[structure][data_dict[structure]['product_coords']])}\n\n")
        for atom, coord in zip(data_dict[structure]['atomsymb'], data_dict[structure][data_dict[structure]['product_coords']]):
            coords = ' '.join([str(float(i)) for i in list(coord)])
            coords = f"{str(atom)} {coords}\n"
            file_3.write(coords)
        file_3.close()
        # Save to rxn_dict.
        rxn_dict[structure] = {'rct_1':name1, 'rct_1_charge':int(np.round(sum(rct1_charge_list))), 'rct_2':name2, 'rct_2_charge':int(np.round(sum(rct2_charge_list))), 'prod':p_name, 'p_charge':data_dict[structure]['charge']}
    
    return rxn_dict

def _xyz_update(structure: str, data_dict: dict, connect_dict: dict, rxn_dict: dict) -> dict:
    """Function to update the .xyz files saved if the initial structure failed.

    Arguments:
    structure (String): A string of the filename that will be used to create the temporary .xyz files for each reactant.
    data_dict (Dictionary): Dictionary with filename (key) and parsed data (value) pairs.
    connect_dict (Dictionary): A dictionary containing the list of atoms that are connected to the reacting atom.
    rxn_dict (Dictionary): A dictionary containing the information on the reactants/products for a given reaction.

    Returns:
    rxn_dict (Dictionary): An updated dictionary containing the information on the reactants/products for a given reaction.

    """
    # Loop through each structure in the connect_dict
    name1 = structure.split('.')[0]+'_reactant_1.xyz'
    file_1 = open(name1, 'w')
    file_1.write(f'{len(connect_dict[structure])}\n\n')
    # Loop through all atoms in the connect_dict for the current structure and pull the coordinates from data_dict
    rct1_charge_list = []
    for atom in connect_dict[structure]:
        coords = list(data_dict[structure][f"max_{data_dict[structure]['reactant_coords']}"][atom])
        coords = ' '.join([str(float(i)) for i in list(coords)])
        rct1_charge_list.append(list(data_dict[structure]['atomcharges'])[atom])
        coords = f"{str(data_dict[structure]['atomsymb'][atom])} {coords}\n"
        file_1.write(coords)
    file_1.close()
    
    # Here we want to get all the atoms in the TS in a list to remove reactant_1 atoms to leave us with reactant_2 atoms
    name2 = structure.split('.')[0]+'_reactant_2.xyz'
    all_atoms = list(range(0, len(data_dict[structure]['coords1'])))
    second_connect = list(set(all_atoms)- set(connect_dict[structure]))

    # Repeat above but this time for reactant_2
    file_2 = open(name2, 'w')
    file_2.write(f'{len(second_connect)}\n\n')
    rct2_charge_list = []
    for atom in second_connect:
        coords = list(data_dict[structure][f"max_{data_dict[structure]['reactant_coords']}"][atom])
        coords = ' '.join([str(float(i)) for i in list(coords)])
        rct2_charge_list.append(list(data_dict[structure]['atomcharges'])[atom])
        coords = f"{str(data_dict[structure]['atomsymb'][atom])} {coords}\n"
        file_2.write(coords)
    file_2.close()

    # Repeat the same for the product
    p_name =  structure.split('.')[0]+'_product.xyz'
    file_3 = open(p_name, 'w')
    file_3.write(f"{len(data_dict[structure][data_dict[structure]['product_coords']])}\n\n")
    for atom, coord in zip(data_dict[structure]['atomsymb'], data_dict[structure][f"max_{data_dict[structure]['product_coords']}"]):
        coords = ' '.join([str(float(i)) for i in list(coord)])
        coords = f"{str(atom)} {coords}\n"
        file_3.write(coords)
    file_3.close()
    # Save to rxn_dict.
    rxn_dict[structure] = {'rct_1':name1, 'rct_1_charge':int(np.round(sum(rct1_charge_list))), 'rct_2':name2, 'rct_2_charge':int(np.round(sum(rct2_charge_list))), 'prod':p_name, 'p_charge':data_dict[structure]['charge']}
    
    return rxn_dict

def _charge_check(rxn_smiles: str) -> bool:
    """Function to check the charges of both reactants and product. Returns False if charge separation has occured in the reaction
    when not wanted.

    Arguments:
    rxn_smiles (String): A string of the reaction SMILES.
    """
    # Split the reaction SMILES back into its components.
    r, p = rxn_smiles.split('>>')[0], rxn_smiles.split('>>')[1]
    r1, r2 = r.split('.')[0], r.split('.')[1]
    # Check for charge in product first
    if '+' in p and '-' in p:
        # Now check if this is charge separation or actually valid.
        if '+' and '-' in r1 or '+' and '-' in r2:
            return True
        else:
            return False
    else:
        return True

def _build_mol(rxn_dict: dict, structure: str) -> Tuple[str, str]:
    """Function to build the reaction SMILES and SMARTS from the temporary .xyz files.

    Arguments:
    rxn_dict (Dictionary): A dictionary containing the information on the reactants/products for a given reaction.
    structure (String): A string of the filename that will be used to create the temporary .xyz files for each reactant.

    Returns:
    rxn_smiles (String): A string containing the reaction SMILES.
    rxn_smarts (String): A string containing the reaction SMARTS.
    """
    # Load in the xyz files.
    r1 = Chem.MolFromXYZFile(rxn_dict[structure]['rct_1'])
    r2 = Chem.MolFromXYZFile(rxn_dict[structure]['rct_2'])
    p = Chem.MolFromXYZFile(rxn_dict[structure]['prod'])

    # Create the mol structures and determine connectivity.
    r1_mol = Chem.Mol(r1)
    rdDetermineBonds.DetermineConnectivity(r1_mol) 
    rdDetermineBonds.DetermineBondOrders(r1_mol, charge=rxn_dict[structure]['rct_1_charge'])   
    r2_mol = Chem.Mol(r2)
    rdDetermineBonds.DetermineConnectivity(r2_mol)    
    rdDetermineBonds.DetermineBondOrders(r2_mol, charge=rxn_dict[structure]['rct_2_charge'])   
    p_mol = Chem.Mol(p)
    rdDetermineBonds.DetermineConnectivity(p_mol)    
    rdDetermineBonds.DetermineBondOrders(p_mol, charge=rxn_dict[structure]['p_charge'])          
    
    # Convert these to SMILES.
    r1_smi = Chem.MolToSmiles(r1_mol)
    r2_smi = Chem.MolToSmiles(r2_mol)
    p_smi = Chem.MolToSmiles(p_mol)

    # Create Reaction SMILES for ReaKE representation. 
    rxn_smiles = f'{r1_smi}.{r2_smi}>>{p_smi}'

    # Convert them to SMARTS
    r1_sma = Chem.MolToSmarts(Chem.MolFromSmiles(r1_smi))
    r2_sma = Chem.MolToSmarts(Chem.MolFromSmiles(r2_smi))
    p_sma = Chem.MolToSmarts(Chem.MolFromSmiles(p_smi))

    # Convert SMARTS to Reaction SMARTS
    rxn_smarts = f'{r1_sma}.{r2_sma}>>{p_sma}'

    return rxn_smiles, rxn_smarts

def _duplicates(filename, seq, rxn_smiles, rxn_smarts) -> bool:
    """Function to check if reaction is already in the file.
    
    Arguments:
    filename (String): A string of the filename being processed.
    seq (String): Either 'w' or 'a' depending on if files already exist.
    rxn_smiles (String): A string containing the reaction SMILES.
    rxn_smarts (String): A string containing the reaction SMARTS.

    Returns:
    True if file/reaction already processed.
    False if a new file/reaction.

    """

    if seq == 'a':
        with open('reaction_smiles.txt') as smi_txt, open('reaction_smarts.txt') as sma_txt:
            if rxn_smiles in smi_txt.read() and rxn_smarts in sma_txt.read():
                print(f'{filename} already processed.')
                return True
            else:
                return False
    else:
        pass
        
def geo(f, k=False, n=False, r=False):
    # Pull out required data from cclib.
    dd = _get_data(f)
    # Determine reacting atoms and get connectivity.
    dd, cd = _get_connected_atoms(dd)
    # Create the temporary xyz files.
    rd = _xyz_creator(dd, cd)

    # Check if files exist and change the sequence to append to existing files.
    seq = ['a' if os.path.isfile('reaction_smiles.txt') and os.path.isfile('reaction_smarts.txt') else 'w'][0]

    # Loop through all available Gaussian TS geometries in rd.
    with open('reaction_smiles.txt', seq) as smi_txt, open('reaction_smarts.txt', seq) as sma_txt:
        for structure in rd.keys():
            # Clean filename.
            filename = os.path.basename(structure)
            try:
                rxn_smiles, rxn_smarts = _build_mol(rd, structure)
                # print('Found Correct TS - 1st Attempt.')
            except:
                try:
                    trd = _xyz_update(structure, dd, cd, rd)
                    rxn_smiles, rxn_smarts = _build_mol(trd, structure)
                    # print(f'{structure} - Found Correct TS - 2nd Attempt.')
                except:
                    rxn_smiles, rxn_smarts = None, None
                    print(f'{filename} - Failed to find suitable transition structure')

            else:
                if rxn_smiles == None and rxn_smarts == None:
                    pass
                else:                
                    # Initial Check that the charges match.
                    if _charge_check(rxn_smiles) == False: # type: ignore
                        print(f"Error in {filename} - split charge in product.")
                    else:
                        # Check if reaction has already been converted previously - save computation.
                        if _duplicates(filename, seq, rxn_smiles, rxn_smarts) is True:
                            pass
                        else:
                        # Write Reaction SMILES and SMARTS to file.
                            if n == False:
                                smi_txt.write(f'{rxn_smiles}\n')
                                sma_txt.write(f'{rxn_smarts}\n')
                            else:
                                smi_txt.write(f'{filename.split('.')[0]},{rxn_smiles}\n')
                                sma_txt.write(f'{filename.split('.')[0]},{rxn_smarts}\n')

                    # Create the reaction.
                    rxn = AllChem.ReactionFromSmarts(rxn_smarts) # type: ignore
                    
                    # Read in any passed arguments.
                    ## Keeping files.
                    if k == False:
                        os.remove(rd[structure]['rct_1'])
                        os.remove(rd[structure]['rct_2'])
                        os.remove(rd[structure]['prod'])
                    ## Saving the reaction as a .png file.
                    if r == True:
                        d2d = Draw.MolDraw2DCairo(800,300)
                        d2d.DrawReaction(rxn)
                        png = d2d.GetDrawingText()
                        open(f'{filename.split(".")[0]}_rxn.png','wb+').write(png)    

    return rxn_smiles, rxn_smarts

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='geomtosmarts', description='Python script for creating Reaction SMILES and SMARTS from optimised TS Gaussian output files.')
    parser.add_argument('-f', dest='filename', required=True)
    parser.add_argument('-k', dest='keep', action='store_true', required=False, help='Argument for keeping the created .xyz files. Defaults as False to not keep them.')
    parser.add_argument('-n', dest='fnames', action='store_true', required=False, help='Argument for saving with filenames. Defaults as False to not save them.')
    parser.add_argument('-r', dest='rxnpng',action='store_true', required=False, help='Argument for saving Reaction SMARTS .png file. Defaults as False to not save the image.')
    args = parser.parse_args()
    geo(args.filename, args.keep, args.fnames, args.rxnpng)