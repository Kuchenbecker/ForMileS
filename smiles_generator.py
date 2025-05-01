import networkx as nx
from rdkit import Chem
from itertools import combinations, product

# Define allowed bonds with maximum bond order
bond_orders = {
    ("C", "C"): [1, 2, 3],  # Single, double, triple bonds
    ("C", "O"): [1, 2],      # Single or double bond
    ("O", "O"): [1],         # Single bond only
}

max_valence = {
    "C": 4,
    "O": 2,
}

def parse_formula(formula):
    """Parses formula like C6O3"""
    import re
    elements = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
    atoms = []
    for elem, count in elements:
        count = int(count) if count else 1
        atoms.extend([elem] * count)
    return atoms

def is_valid(graph):
    """Check if all atoms have valid valence"""
    for node in graph.nodes:
        atom = graph.nodes[node]['element']
        valence = sum(data['order'] for _, _, data in graph.edges(node, data=True))
        if valence > max_valence[atom]:
            return False
    return True

def atoms_graph(atoms):
    """Create a graph with atoms as nodes"""
    G = nx.Graph()
    for i, atom in enumerate(atoms):
        G.add_node(i, element=atom)
    return G

def generate_graphs(atoms):
    """Generate all possible graphs with bond multiplicities, allowing cycles"""
    G = atoms_graph(atoms)
    n = len(atoms)
    graphs = []

    possible_edges = list(combinations(range(n), 2))

    # Allow more than tree edges to enable cycles
    for num_edges in range(n - 1, len(possible_edges) + 1):
        for selected_edges in combinations(possible_edges, num_edges):
            g = G.copy()
            valid = True
            for (i, j) in selected_edges:
                a1, a2 = g.nodes[i]['element'], g.nodes[j]['element']
                allowed_orders = bond_orders.get((a1, a2)) or bond_orders.get((a2, a1))
                if not allowed_orders:
                    valid = False
                    break
                g.add_edge(i, j, order=1)
            if valid and nx.is_connected(g):
                expanded_graphs = expand_bond_orders(g)
                for eg in expanded_graphs:
                    if is_valid(eg):
                        graphs.append(eg)
    return graphs

def expand_bond_orders(graph):
    """Expand graphs by trying higher bond orders"""
    graphs = []
    edges = list(graph.edges(data=True))

    bond_options = []
    for (i, j, data) in edges:
        a1 = graph.nodes[i]['element']
        a2 = graph.nodes[j]['element']
        allowed_orders = bond_orders.get((a1, a2)) or bond_orders.get((a2, a1))
        bond_options.append([(i, j, order) for order in allowed_orders])

    for bond_combination in product(*bond_options):
        g = graph.copy()
        for (i, j, order) in bond_combination:
            g[i][j]['order'] = order
        graphs.append(g)
    return graphs

def graph_to_rdkit_mol(graph):
    """Convert NetworkX graph to RDKit Mol (without hydrogens)"""
    rw_mol = Chem.RWMol()
    node_to_idx = {}

    for node in graph.nodes:
        atom = Chem.Atom(graph.nodes[node]['element'])
        idx = rw_mol.AddAtom(atom)
        node_to_idx[node] = idx

    for i, j, data in graph.edges(data=True):
        rw_mol.AddBond(node_to_idx[i], node_to_idx[j], number_to_bondtype(data['order']))

    return rw_mol

def number_to_bondtype(order):
    if order == 1:
        return Chem.BondType.SINGLE
    elif order == 2:
        return Chem.BondType.DOUBLE
    elif order == 3:
        return Chem.BondType.TRIPLE
    else:
        raise ValueError("Unsupported bond order")

def generate_smiles(formula):
    atoms = parse_formula(formula)
    graphs = generate_graphs(atoms)
    smiles_list = []
    for g in graphs:
        mol = graph_to_rdkit_mol(g)
        try:
            smiles = Chem.MolToSmiles(mol, canonical=True)
            smiles_list.append(smiles)
        except:
            continue
    return list(set(smiles_list))

if __name__ == "__main__":
    formula = "C3O"
    smiles_list = generate_smiles(formula)
    print("Generated SMILES:")
    for s in smiles_list:
        print(s)

