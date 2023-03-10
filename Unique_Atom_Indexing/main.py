from pymol import cmd
from rdkit import Chem
import morgan


def main():
    '''
    Run everything
    :return: none
    '''

    input = '../brom_bb_360_mp2.sdf'
    molecules = []
    with Chem.SDMolSupplier(input, sanitize = False, removeHs=False) as supply:
        molecules = [mol for mol in supply]
    test = morgan.Molecule('test',molecules[1])
    print(list(molecules[1].GetPropNames()))
    test.create_graph()
    print(test.graph)
    test.morgan()
    test.draw_graph('unique_index')
    #test.draw_graph()
    #### test differnt stuff








main()