import networkx as nx
import numpy as np
import rule_set
import matplotlib.pyplot as plt
from rdkit import Chem
from functools import cmp_to_key


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
        self.ambiguity_rules = rule_set.rule_dict()

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

    def relax(self, max_cycles =50000):
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

        # loop unit less unique lables or max cycles
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

    def resolve_ambiguities(self, neighbors, current_node) -> list :
        '''
        Sort the list of nodes according to the rules provided.
        Applies all rules starting at the first dictionary entry.
        Additionally, a distance based function will be applied as rule
        :param node: current node
        :param neighbors: current nodes neighbors as list
        :return: sorted list of nodes
        '''
        # generate attributes for each node
        # stable sort by given rules
        # stable sort by distance
        # return list
        rules = list(self.ambiguity_rules.values())
        G = self.graph.nodes

        # sort neighbors by last rule first
        for rule in rules[::-1]:
            def rule_wrapper(x, y):
                # wrapper to pars current_node argument  into cmp_to_key function
                return rule(self, node1=x, node2=y, V=current_node, center = (0,0,0))
            neighbors.sort(key = cmp_to_key(rule_wrapper), reverse=True)

        # ensure sorting by EC lable
        neighbors.sort(key=lambda x: G[x]['EC'], reverse=True)
        return neighbors


    def subgraph_morgan(self, subset, C_start=0):
        '''
        Performs Morgan algo for canonical enumeration on a subset of the graph
        An alternative set of rules to solve ambiguities can be passed as argument.
        In case ambiguities can not be solved, enumeration is solved via a distance approach.
        :param subset: list of nodes describing the subgraph
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
        for node in G.nodes:
            G.nodes[node]['unique_index'] = C_start-1
        EC_lables = dict(G.nodes(data='EC'))

        # get node with highest EC
        V = max(EC_lables, key = EC_lables.get)
        V_queue = []
        G.nodes[V]['unique_index'] = C_start
        C = C_start+1

        # repeat until all index are not initialisation value
        while C_start-1 in dict(G.nodes(data='unique_index')).values():
            neighbors = [n for n in list(nx.neighbors(G,V)) if (G.nodes[n]['unique_index'] == C_start-1)]
            neighbors.sort(key=lambda x: G.nodes[x]['EC'], reverse=True)

            # resolve ambiguities
            if len(neighbors) != len(set([G.nodes[n]['EC'] for n in neighbors])):
                neighbors = self.resolve_ambiguities(neighbors, V)

            # assige unique index
            for n in neighbors:
                G.nodes[n]['unique_index'] = C
                V_queue.append(n)
                C += 1

            # next node
            if len(V_queue) == 0:
                break
            V = V_queue.pop(0)


    def morgan(self):
        '''
        Performs Morgan algo for canonical enumeration on a graph
        An alternative set of rules to solve ambiguities can be passed as argument.
        In case ambiguities can not be solved, enumeration is solved via a distance approach.
        :param rules: dictionary containing an alternative set of rules for ambiguity resolution.
        :return: none
        '''
        # run all morgan algo subsections
        #   run relax on full graph
        #   find subgraphs
        #       run subgrap_morgan on each subgaph
        #   return mapping of atom index to unique_index
        self.relax()
        nodes = self.graph.nodes
        start = 1
        for g in sorted(nx.connected_components(self.graph), key = len, reverse=True):
            print(g)
            self.subgraph_morgan(g,C_start= start)
            #start = start + len(g)

        for node in nodes:
            self.mapping[nodes[node]['original_atom_idx']] = nodes[node]['unique_index']
        print(self.mapping)






