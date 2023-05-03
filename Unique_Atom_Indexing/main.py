from rdkit import Chem
import morgan

def mapping_sdf(file_path: str, restart_index: (True, False) = True, method: str = 'morgan') -> dict:
    '''
    Atom index mapping for all molecules in sdf style files. Returns the mapping of the molecule as dictionary of
    lists. A list is keyed to the molecule name with atom i being mapped to entry A[i].
    :param file_path: str, path to sdf molecule file
    :param method: method used to create mapping, one of 'morgan, pairs or spanning_tree'
    :param restart_index: bool, restart the atom mapping at 0 if there are multiple molecules in molecule block.
    Otherwise, continue mapping molecules in object with consecutive indices. Mapping ist not consecutive across
    molecule blocks in sdf file.
    :return: dictionary of lists, dic{'mol name': A[i], ...},
    '''
    molecules = {}
    # read file
    with Chem.ForwardSDMolSupplier(file_path, sanitize = False, removeHs=False) as supply:
        # loop mol blocks
        for mol in supply:
            if mol is None: continue
            # create object and use enumeration method
            name = mol.GetProp('_Name')
            molecule_object = morgan.Molecule.from_Mol(name,mol)
            if method == 'morgan':
                molecule_object.morgan(reset=restart_index)
            elif method == 'pairs':
                molecule_object.pairs_method(reset=restart_index)
            elif method == 'spanning_tree':
                molecule_object.spanning_tree_method(reset=restart_index)
            else:
                raise Exception(f'Method {method} is not a viable method')
            molecules[name] = molecule_object.get_mapping()
    return molecules

def mapping_xyz(file_path: str, restart_index: (True, False) = True, method: str = 'morgan') -> dict:
    '''
    Atom index mapping for all molecules in xyz style files. Returns the mapping of the molecule as dictionary of
    lists. A list is keyed to the molecule name with atom i being mapped to entry A[i].
    :param file: string, path to file
    :param restart_index: bool, wether to restart index for multiple molecules in a single mol block.
    :param method: str, Enumeration method to use, one of morgan, pairs or spanning_tree.
    :return: dictionary of lists, dic{'mol name': A[i], ...},
    '''
    xyz_cords = []
    number_atoms = 0
    mol_name = ''
    molecules = {}
    with open(file_path, 'r') as xyz_file:
        # read file line wise
        for line_number, line in enumerate(xyz_file):
            if line_number == 0:
                numer_atoms = int(line)
            elif line_number == 1:
                mol_name = line.split()[0]
            else:
                element, x, y, z, = line.split()
                xyz_cords.append((element, float(x), float(y), float(z)))
        # create molecule object and run enumeration
        molecule_object = morgan.Molecule.from_XYZ(mol_name, xyz_cords)
        if method == 'morgan':
            molecule_object.morgan(reset=restart_index)
        elif method == 'pairs':
            molecule_object.pairs_method(reset=restart_index)
        elif method == 'spanning_tree':
            molecule_object.spanning_tree_method(reset=restart_index)
        else:
            raise Exception(f'Method {method} is not a viable method')
        molecules[mol_name] = molecule_object.get_mapping()
    return molecules

    ####teststts

print(mapping_xyz('glucose.xyz', method='pairs'))