import networkx as nx
import math
import os
from itertools import combinations, product
from rdkit import Chem

bond_orders = {
    ("C", "C"): [1, 2, 3],
    ("C", "O"): [1, 2],
    ("O", "O"): [1],
    ("C", "N"): [1, 2, 3],
    ("N", "N"): [1, 2, 3],
}

max_valence = {
    "C": 4,
    "O": 2,
    "N": 3,
}

def parse_formula(FORMULA):
    import re
    elements = re.findall(r'([A-Z][a-z]*)(\d*)', FORMULA)
    atoms = []
    for elem, count in elements:
        count = int(count) if count else 1
        atoms.extend([elem] * count)
    return atoms

def is_valid(graph):
    for node in graph.nodes:
        atom = graph.nodes[node]['element']
        valence = sum(data['order'] for _, _, data in graph.edges(node, data=True))
        if valence > max_valence[atom]:
            return False
    return True

def atoms_graph(atoms):
    G = nx.Graph()
    for i, atom in enumerate(atoms):
        G.add_node(i, element=atom)
    return G

def expand_bond_orders(graph):
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

def number_to_bondtype(order):
    if order == 1: return Chem.BondType.SINGLE
    if order == 2: return Chem.BondType.DOUBLE
    if order == 3: return Chem.BondType.TRIPLE
    raise ValueError("Invalid bond order")

def graph_to_rdkit_mol(graph):
    rw_mol = Chem.RWMol()
    node_to_idx = {}
    for node in graph.nodes:
        atom = Chem.Atom(graph.nodes[node]['element'])
        idx = rw_mol.AddAtom(atom)
        node_to_idx[node] = idx
    for i, j, data in graph.edges(data=True):
        rw_mol.AddBond(node_to_idx[i], node_to_idx[j], number_to_bondtype(data['order']))
    return rw_mol

def generate_graphs_lazy(atoms, should_stop=lambda: False, update_callback=None):
    G_base = atoms_graph(atoms)
    n = len(atoms)
    possible_edges = list(combinations(range(n), 2))
    total_combinations = math.comb(len(possible_edges), n - 1)
    processed = 0
    last_update = 0

    for edge_comb in combinations(possible_edges, n - 1):
        if should_stop():
            break
            
        g = G_base.copy()
        valid = True

        for (i, j) in edge_comb:
            a1, a2 = g.nodes[i]['element'], g.nodes[j]['element']
            if (a1, a2) in bond_orders or (a2, a1) in bond_orders:
                g.add_edge(i, j, order=1)
            else:
                valid = False
                break

        if valid and nx.is_tree(g) and is_valid(g):
            expanded = expand_bond_orders(g)
            for eg in expanded:
                if is_valid(eg):
                    yield eg
                    
        processed += 1
        if update_callback and (processed % 10 == 0 or processed == total_combinations):
            update_callback(processed, total_combinations)

def generate_smiles(FORMULA, CHARGE, output_file=None, should_stop=lambda: False, update_callback=None):
    output_dir = f"OutputFiles_{FORMULA}_Charge_{CHARGE}"
    output_path = os.path.join(output_dir, output_file or f"nSMILES_{FORMULA}.txt")
    atoms = parse_formula(FORMULA)
    graphs = list(generate_graphs_lazy(atoms, should_stop, update_callback))
    total = len(graphs)
    smiles_set = set()

    for i, g in enumerate(graphs):
        if should_stop():
            print("Stopped during SMILES generation.")
            break
        mol = graph_to_rdkit_mol(g)
        try:
            smiles = Chem.MolToSmiles(mol, canonical=True)
            smiles_set.add(smiles)
        except:
            continue
        if update_callback and (i % 10 == 0 or i == total - 1):
            update_callback(i + 1, total)

    with open(output_path, "w") as f:
        for smiles in sorted(smiles_set):
            f.write(smiles + "\n")
    print(f"Saved {len(smiles_set)} SMILES to '{output_path}'")