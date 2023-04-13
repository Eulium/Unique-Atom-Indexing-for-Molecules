import itertools
import rule_set
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from rdkit import Chem
from functools import cmp_to_key
from scipy.spatial.distance import cdist
from scipy.stats import rankdata

class Molecule:
    def __init__(self, mol_name: str, molecule_object, file_type: str):
        '''
        Initialises as Molecule class
        :param molecule_object: a RDkit molecule object
        :param mol_name: name of the molecule
        '''
        self.original_molecule = molecule_object
        self.name = mol_name
        self.file_type = file_type
        self.ambiguity_rules = rule_set.rule_dict()
        self.graph = nx.Graph()
        self.d_matrix = None
        if file_type == 'Mol':
            self.graph = self.create_graph_from_mol()
            self.mapping = [None] * len(molecule_object.GetAtoms())
        elif file_type == 'xyz':
            self.graph = self.create_graph_from_xyz()
            self.mapping = [None] * len(molecule_object)
    @classmethod
    def from_Mol(cls, mol_name: str,  rdchem_mol_object: Chem.rdchem.Mol):
        molecule = cls(mol_name, rdchem_mol_object, 'Mol')
        return molecule

    @classmethod
    def from_XYZ(cls, mol_name: str, xyz_coords: list):
        molecule =cls(mol_name, xyz_coords, 'xyz')
        return molecule

    def create_graph_from_mol(self):
        '''
        Create an undirected nx graph from a Molecule object
        Nodes contain atom object and original atom index
        :return: none
        '''
        mol = self.original_molecule
        graph = nx.Graph()
        atoms = list(mol.GetAtoms())
        nodes = [{'original_atom_idx': atom.GetIdx(),
                  'element': atom.GetSymbol(),
                  'pos': mol.GetConformer().GetAtomPosition(atom.GetIdx()),
                  'atom_object': atom} for atom in atoms]
        edges = [(atom.GetIdx(), ng.GetIdx(),
                  {'bond_type': str(mol.GetBondBetweenAtoms(atom.GetIdx(), ng.GetIdx()).GetBondType())})
                 for atom in atoms
                 for ng in atom.GetNeighbors()
                 ]
        nodes_to_add = [(node['original_atom_idx'], node) for node in nodes]
        graph.add_nodes_from(nodes_to_add)
        graph.add_edges_from(edges)
        return graph

    def create_graph_from_xyz(self):

        g = nx.Graph()
        mol = self.original_molecule
        coords = [(xyz[1:]) for xyz in mol]
        dist_matrix = cdist(np.array(coords), np.array(coords))
        self.d_matrix = dist_matrix
        nodes = [{'original_atom_idx': num,
                  'element': atom,
                  'pos': (x,y,z),
                  'atom_object': (atom, x, y, z)}
                 for num,(atom,x,y,z) in enumerate(mol)]
        edges = []
        for x in range(0, len(mol)):
            for y in range(x+1, len(mol)):
                if dist_matrix[x][y] <= 2:
                    edges.append((x, y, {'bond_type':'SINGLE'}))

        g.add_nodes_from([(node['original_atom_idx'], node) for node in nodes])
        g.add_edges_from(edges)
        return g
    def get_mapping(self):
        '''
        Returns the mapping of the original atom enumeration to unique labeling.
        A[i] = unique label of original atom i
        :return: [int, ...] list of int,
        '''
        if self.mapping[0] == None:
            g = self.graph
            if not 'unique_index' in g.nodes[-1].keys():
                raise Exception('No index created yet, pls run morgan, pairs_method or other first')
            else:
                self.update_mappin()
                return self.mapping
        else:
            return self.mapping
    def update_mappin(self):
        '''
        Updates the molecules mapping attribute.
        :return:
        '''
        g = self.graph
        mapping = [g.nodes[i]['unique_index'] for i in g.nodes]
        self.mapping = mapping

    def draw_graph(self, label_key = False):
        '''
        Draw and show the molecule graph
        :param label_key: key for label dictionary
        :return:
        '''
        if label_key != False:
            nx.draw_networkx(self.graph, labels=dict(self.graph.nodes(data=label_key, default="Not Available")))
        else:
            nx.draw_networkx(self.graph)
        plt.show()


    def center_pos(self, subset: nx.Graph):
        '''
        For a given graph subset calculate the center via the atoms mean position
        :param subset: networkx subgraph
        :return: (float, float, float), center as (x, y, z) triple
        '''
        center = [(self.graph.nodes[n]['pos']) for n in subset]
        center = [sum(center[i]) / len(center[i]) for i in range(0, 3)]
        return center

    def relax(self, max_cycles =50000):
        '''
        Does the Morgan relaxation step for a graph instance.
        :param max_cycles: max amount of relaxation cylces for latest termination
        :return: none
        '''
        g = self.graph
        for node in g.nodes:
            g.nodes[node]['EC_l'] = 1
            g.nodes[node]['EC_c'] = 1
        number_l_labels = 1
        number_c_labels = 1
        cycle = 0

        # loop unit less unique labels or max cycles
        while cycle <= max_cycles:
            for node in g.nodes:
                g.nodes[node]['EC_c'] = sum([g.nodes(data= 'EC_l')[n] for n in g.neighbors(node)])

            number_c_labels = len(set(dict((g.nodes(data='EC_c'))).values()))
            # set current EC labels as last rounds labels
            if number_c_labels <= number_l_labels:
                for node in g.nodes:
                    g.nodes[node]['EC_c'] = g.nodes[node]['EC_l']
                break
            # replace last rounds label with current round labels
            else:
                for node in g.nodes:
                    g.nodes[node]['EC_l'] = g.nodes[node]['EC_c']
                number_l_labels = number_c_labels
            cycle += 1

        # remove unnecessary attributes
        for node in g.nodes:
            g.nodes[node]['EC'] = g.nodes[node]['EC_c']
            del g.nodes[node]['EC_l']
            del g.nodes[node]['EC_c']

    def resolve_ambiguities(self, neighbors, current_node, center) -> list :
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
                return rule(self, node1=x, node2=y, V=current_node, center = center)
            neighbors.sort(key = cmp_to_key(rule_wrapper), reverse=True)

        # ensure sorting by EC label
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
        center = self.center_pos(subset)
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
                neighbors = self.resolve_ambiguities(neighbors, V, center)

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
        start = 0
        for g in sorted(nx.connected_components(self.graph), key = len, reverse=True):
            self.subgraph_morgan(g,C_start= start)
            start = start + len(g)
        self.update_mappin()


    def subgraph_pairs_methods(self, subset, c_start: (int, float)):
        '''
        Uses the distance between pairs of atoms to create a unique numbering.
        Algo:
            Given a set of nodes, with xyz coords
            1. fully connect nodes
            2. label edges with distance between points.
            3. sort edges by weight, low to high
            4. assign an ordering to edges based on edge weight
            5. label nodes as a string of their edge labels
            6. sort node by labels,
            7. assign rank as new label
        :param subset: subset of nodes to apply algo to, node i is also atom i in xyz file
        :return : list, [(node, unique_label),...], return a list of tuples with node numer and unique label number.
        '''
        # get distance matrix and create a graph subset
        d_matrix = self.d_matrix
        g = nx.Graph()
        g.add_nodes_from(subset)
        # get distances for graph subset only, sort by distance
        weights = [(x, y, d_matrix[x][y]) for x, y in itertools.combinations(subset, 2)]
        weights.sort(key = lambda x: x[2], reverse=False)
        # loop all determined weights, assign weights and rank to edegs
        for rank, (x, y, dist) in enumerate(weights):
            # edge no present
            if not (x, y) in g.edges or not (y, x) in g.edges:
                g.add_edge(x, y, weight = dist)
                g.edges[x, y]['rank'] = rank
            # edge present
            else:
                g.edges[x, y]['weight'] = dist
                g.edges[x, y]['rank'] = rank
        # for every node create new label as sorted list of edge ranks
        for node_a in g.nodes:
            node_a_edges = [(x, y) for x, y, dist in weights if (x == node_a) or (y == node_a)]
            labels = [g.edges[x, y]['rank'] for x, y in node_a_edges]
            g.nodes[node_a]['rank_label'] = labels
        # sort nodes by rank labels, translate into new unique label
        labels = [(x, g.nodes[x]['rank_label']) for x in g.nodes]
        labels.sort(key=lambda x: x[1])

        # list of nodes and their new label
        re = [(node, rank+c_start) for rank, (node, label_list) in enumerate(labels)]
        a = rankdata(labels, method='min', axis=0)
        print(a)
        return re

    def pairs_method(self):
        '''
        Create unique atom enumeration for the Molecule object. Uses as pairwise distances method instead of Morgan algo.
        In the case of multiple disconnected components in the Molecule graph creates unique labels for all of them.
        Works for molecules generated from Mol and XYZ file types.
        :return: none
        '''
        start = 0
        for subset in sorted(nx.connected_components(self.graph), key = len, reverse=True):
            mapping_subset = self.subgraph_pairs_methods(subset, c_start= start)
            for atom_number, atom_mapping in mapping_subset:
                self.graph.nodes[atom_number]['unique_index'] = atom_mapping
            start = start + len(subset)





