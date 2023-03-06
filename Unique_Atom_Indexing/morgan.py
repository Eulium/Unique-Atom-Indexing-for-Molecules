import networkx as nx
import matplotlib.pyplot as plt
from rdkit import Chem


class Molecule:


    def __init__(self, mol_name: str, molecule_object: Chem.rdchem.Mol):
        '''
        Initialises as Molecule class
        :param molecule_object: a RDkit molecule object
        :param mol_name: name of the molecule
        '''
        self.original_molecule = molecule_object
        self.graph = nx.Graph()
        self.name = mol_name
        self.mapping = [0]*len(molecule_object.GetAtoms())
        self.ambiguity_rules = {'bond_type': ['triple', 'double', 'single'],
                           'element_type': ['O', 'C' 'S', 'N', 'H']}

    def draw_graph(self, label_key = False):
        '''
        Draw and show the molecule graph
        :param label_key: key for label dictionary
        :return:
        '''
        if bool(label_key):
            nx.draw_networkx(self.graph, labels=dict(self.graph.nodes(data=label_key, default="Not Available")))
        else:
            nx.draw_networkx(self.graph)
        plt.show()


    def create_graph(self):
        '''
        Create an undirected nx graph from a Molecule object
        Nodes contain atom object and original atom index
        :return: none
        '''
        mol = self.original_molecule
        graph = self.graph
        atoms = list(mol.GetAtoms())
        nodes = [{'original_atom_idx': atom.GetIdx(), 'atom_object': atom} for atom in atoms]
        edges = [ (atom.GetIdx(), ng.GetIdx(), {'bond_type':str(mol.GetBondBetweenAtoms(atom.GetIdx(), ng.GetIdx()).GetBondType())})
                 for atom in atoms
                 for ng in atom.GetNeighbors()
                  ]
        graph.add_nodes_from([(node['original_atom_idx'], node) for node in nodes])
        graph.add_edges_from(edges)

    def relax(self, max_cycles =5000):
        '''
        Does the Morgan relaxation step for a graph instance.
        :param max_cycles: max amount of relaxation cylces for latest termination
        :return: none
        '''
        G = self.graph
        for node in G.nodes:
            G.nodes[node]['EC_l'] = 1
            G.nodes[node]['EC_c'] = 1
        number_L_labels = 1
        number_C_labels = 1
        cycle = 0

        while cycle <= max_cycles:
            for node in G.nodes:
                G.nodes[node]['EC_c'] = sum([G.nodes(data= 'EC_l')[n] for n in G.neighbors(node)])

            number_C_labels = len(set(dict((G.nodes(data='EC_c'))).values()))
            # set current EC labels as last rounds labels
            if number_C_labels <= number_L_labels:
                for node in G.nodes:
                    G.nodes[node]['EC_c'] = G.nodes[node]['EC_l']
                break
            # replace last rounds label with current round labels
            else:
                for node in G.nodes:
                    G.nodes[node]['EC_l'] = G.nodes[node]['EC_c']
                number_L_labels = number_C_labels
            cycle += 1

        # remove unnecessary attributes
        for node in G.nodes:
            G.nodes[node]['EC'] = G.nodes[node]['EC_c']
            del G.nodes[node]['EC_l']
            del G.nodes[node]['EC_c']

    def resolve_ambiguities(self, nodes, rules):
        '''
        Sort the list of nodes according to the rules provided.
        Applies all rules starting at the first dictionary entry.
        Additionaly a distance based function will be applied as rule
        :param nodes: list of nodes
        :param rules: a dict of rules
        :return: sorted list of nodes
        '''


    def subgraph_morgan(self, subset, rules=False, C_start=0):
        '''
        Performs Morgan algo for canonical enumeration on a subset of the graph
        An alternative set of rules to solve ambiguities can be passed as argument.
        In case ambiguities can not be solved, enumeration is solved via a distance approach.
        :param subset: list of nodes describing the subgraph
        :param rules: dictionary containing an alternative set of rules for ambiguity resolution.
        :param C_start: int which denotes the lowest enumeration
        :return: none
        '''
        # Morgan algo:
        # initiate unique index 0
        # set node with max EC (V) to unique index 1
        # until all unique index are not 0
        #   get all neighbors of V:
        #       if unique EC labels: enumerate
        #       else: resolve ambiguity and enumerate

        # initiate unique index
        G = self.graph.subgraph(subset)
        print(G)
        print(G.nodes)
        for node in G.nodes:
            G.nodes[node]['unique_index'] = C_start-1
        EC_lables = dict(G.nodes(data='EC'))

        # get node with highest EC
        V = max(EC_lables, key = EC_lables.get)
        V_queue = []
        G.nodes[V]['unique_index'] = C_start
        C = C_start+1

        # repeat until all index are not 0
        while C_start-1 in dict(G.nodes(data='unique_index')).values():
            neighbors = [n for n in list(nx.neighbors(G,V)) if (G.nodes[n]['unique_index'] == C_start-1)]
            neighbors.sort(key=lambda x: G.nodes[x]['EC'], reverse=True)

            if len(neighbors) != len(set([G.nodes[n]['EC'] for n in neighbors])):
                # resolve ambiguities
                for n in neighbors:
                    G.nodes[n]['unique_index'] = C
                    V_queue.append(n)
                    C += 1
            else:
                # no ambiguities
                for n in neighbors:
                    G.nodes[n]['unique_index'] = C
                    V_queue.append(n)
                    C += 1
            #print(f'node {V}, neighbors:{neighbors}, queue {V_queue}')
            # next node
            if len(V_queue) == 0:
                break
            V = V_queue.pop(0)


    def morgan(self, rules = False):
        '''
        Performs Morgan algo for canonical enumeration on a graph
        An alternative set of rules to solve ambiguities can be passed as argument.
        In case ambiguities can not be solved, enumeration is solved via a distance approach.
        :param rules: dictionary containing an alternative set of rules for ambiguity resolution.
        :return: none
        '''

        self.relax()
        self.subgraph_morgan(range(0,7), self.ambiguity_rules)
        #self.subgraph_morgan(range(7,12))






