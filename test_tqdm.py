from rdkit import Chem
from rdkit.Chem import rdmolops
import networkx as nx
from itertools import combinations
import re
from tqdm import tqdm
import multiprocessing as mp

valences = {'C': 4, 'H': 1, 'O': 2}
bond_orders = {('C', 'C'): [1, 2, 3],
               ('C', 'H'): [1],
               ('C', 'O'): [1, 2],
               ('O', 'H'): [1],
               ('O', 'O'): [1, 2]}

def parse_formula(formula):
    tokens = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
    atoms = []
    for elem, count in tokens:
        count = int(count) if count else 1
        atoms.extend([elem] * count)
    return atoms

def atoms_graph(atoms):
    G = nx.Graph()
    for i, atom in enumerate(atoms):
        G.add_node(i, element=atom)
    return G

def contains_cycle(graph):
    try:
        nx.find_cycle(graph)
        return True
    except nx.NetworkXNoCycle:
        return False

def generate_graph_candidates(atoms):
    G = atoms_graph(atoms)
    n = len(atoms)
    candidates = []
    possible_edges = list(combinations(range(n), 2))

    for selected_edges in combinations(possible_edges, n - 1):
        g = G.copy()
        valid = True
        for (i, j) in selected_edges:
            a1, a2 = g.nodes[i]['element'], g.nodes[j]['element']
            allowed_orders = bond_orders.get((a1, a2)) or bond_orders.get((a2, a1))
            if not allowed_orders:
                valid = False
                break
            g.add_edge(i, j, order=1)
        if valid and nx.is_connected(g) and not contains_cycle(g):
            candidates.append(g)
    return candidates

def expand_bond_orders(G):
    graphs = [G.copy()]
    for (i, j) in G.edges:
        a1, a2 = G.nodes[i]['element'], G.nodes[j]['element']
        allowed_orders = bond_orders.get((a1, a2)) or bond_orders.get((a2, a1))
        new_graphs = []
        for g in graphs:
            for order in allowed_orders:
                g2 = g.copy()
                g2[i][j]['order'] = order
                new_graphs.append(g2)
        graphs = new_graphs
    return graphs

def is_valid(G):
    valence = {i: 0 for i in G.nodes}
    for i, j in G.edges:
        order = G[i][j]['order']
        valence[i] += order
        valence[j] += order
    for i in G.nodes:
        atom = G.nodes[i]['element']
        if valence[i] > valences[atom]:
            return False
    return True

def graph_to_rdkit_mol(G):
    mol = Chem.RWMol()
    idx_map = {}
    for i in G.nodes:
        a = Chem.Atom(G.nodes[i]['element'])
        idx = mol.AddAtom(a)
        idx_map[i] = idx
    for i, j in G.edges:
        order = G[i][j]['order']
        mol.AddBond(idx_map[i], idx_map[j], Chem.BondType.values[order])
    rdmolops.SanitizeMol(mol)
    return mol

def expand_and_filter_graph(g):
    for eg in expand_bond_orders(g):
        if is_valid(eg):
            try:
                mol = graph_to_rdkit_mol(eg)
                smi = Chem.MolToSmiles(mol, canonical=True)
                return smi
            except:
                return None
    return None

def generate_smiles(formula):
    atoms = parse_formula(formula)
    graph_candidates = generate_graph_candidates(atoms)

    with mp.Pool(processes=mp.cpu_count()) as pool:
        results = list(tqdm(pool.imap(expand_and_filter_graph, graph_candidates),
                            total=len(graph_candidates),
                            desc="Generating SMILES"))

    smiles_list = list(set(filter(None, results)))
    return smiles_list

if __name__ == "__main__":
    formula = "C2H6O"
    smiles = generate_smiles(formula)
    print(f"Generated {len(smiles)} unique SMILES:")
    for s in smiles:
        print(s)

