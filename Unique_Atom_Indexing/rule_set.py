'''
 Use as interface to create ambiguity rules as a comparison function.
 Right functions here and add 'string:func_name' to rule_dict
 Each rule will receive 3 key word arguments, node1, node 2 and the current_node.
 Node1 and Node2 are the two atoms to be compared, while current_node is the morgan algo current node for enumeration.
 A comparison function must return -1 if node1 is to be order lower than node2, 1 if node2 < node1
 and 0 if they are of the same order.
 self.graph.nodes
'''

def rule_dict():
    # returns ambiguity rules as a dictionary
    rules = {'distance': by_distance, 'bond_type': bond_type, 'elm_type': elm_type}
    return rules

def elm_type(self, **kwargs):
    # elemet type wise comparison
    # node1 = x, node2 = y, V = current_node
    Nodes = self.graph.nodes
    Edges = self.graph.edges
    ord = ['H', 'B', 'C', 'O', 'N', 'S', 'Cl']
    a = Nodes[kwargs['node1']]['atom_object'].GetSymbol()
    b = Nodes[kwargs['node2']]['atom_object'].GetSymbol()
    a = len(ord) if a not in ord else ord.index(a)
    b = len(ord) if b not in ord else ord.index(b)
    re = -1 if a < b else 1 if b < a else 0
    return re

def bond_type(self, **kwargs):
    # bond type wise comparison
    Edges = self.graph.edges
    order = ['SINGLE', 'DOUBLE', 'TRIPLE', 'QUADRUPLE', 'QUINTUPLE', 'HEXTUPLE', 'AROMATIC', 'IONIC']
    bond_a = Edges[kwargs['node1'], kwargs['V']]['bond_type']
    bond_b = Edges[kwargs['node2'], kwargs['V']]['bond_type']
    a = len(order) if not bond_a in order else order.index(bond_a)
    b = len(order) if not bond_b in order else order.index(bond_b)
    re = -1 if a < b else 1 if b < a else 0
    return re

def by_distance(self, **kwargs):
    # compares nodes based on the distance to a center
    Nodes = self.graph.nodes
    center = kwargs['center']
    pos_a = self.original_molecule.GetConformer().GetAtomPosition(Nodes[kwargs['node1']]['original_atom_idx'])
    pos_b = self.original_molecule.GetConformer().GetAtomPosition(Nodes[kwargs['node2']]['original_atom_idx'])
    def dist(a, b):
        d = sum([(x-y)**2 for x, y in zip(a,b)])**(0.5)
        return d
    d1 = dist(pos_a, center)
    d2 = dist(pos_b, center)
    re = -1 if d1 < d2 else 1 if d2 < d1 else 0
    return re