from pymol import cmd
from rdkit import Chem
import morgan


def main_sdf():
    '''
    Run everything
    :return: none
    '''
    input = '../brom_bb_360_mp2.sdf'
    molecules = []
    with Chem.SDMolSupplier(input, sanitize = False, removeHs=False) as supply:
        molecules = [mol for mol in supply]
    brom_pair = morgan.Molecule.from_Mol('brom_pair',molecules[1])


    # run morgan algo
    brom_pair.morgan()
    brom_pair.draw_graph('unique_index')

def read_xyz(file_path):
    '''
    Read a molecule's xyz file.
    :param file: string, path to file
    :return: string, [(string, float, float, float), ...] molecule name and list of tuples with element name and xyz coords
    '''
    xyz_cords = []
    with open(file_path, 'r') as xyz_file:
        for line_number, line in enumerate(xyz_file):
            if line_number == 0:
                numer_atoms = int(line)
            elif line_number == 1:
                mol_name = line.split()[0]
            else:
                element, x, y, z, = line.split()
                xyz_cords.append((element, float(x), float(y), float(z)))
    return mol_name, xyz_cords


mol_name, xyz = read_xyz('glucose.xyz')
xyz = morgan.Molecule.from_XYZ(mol_name, xyz)
xyz.draw_graph(label_key='original_atom_idx')
#xyz.morgan()

xyz.pairs_method()
xyz.draw_graph(label_key='unique_index')

#main_sdf()