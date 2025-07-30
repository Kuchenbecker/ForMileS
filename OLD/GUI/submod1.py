#########################################################################
# This module generates molecular graphs using graph theory and          #
# NetworkX, with multiprocessing support for improved performance       #
#########################################################################

import networkx as nx
import math
import os
from functools import lru_cache
from multiprocessing import Pool, cpu_count
from itertools import combinations, product

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from tqdm import tqdm 

############ Bond and Valence Rules ###########
bond_orders = {
    ("C", "C"): [1, 2],
    ("C", "O"): [1, 2],
    ("O", "O"): [1],
}

max_valence = {
    "C": 4,
    "O": 2,
}

####### Functions for SMILES Generation #######

def parse_formula(FORMULA):
    """Parse chemical formula into list of atoms"""
    import re
    elements = re.findall(r'([A-Z][a-z]*)(\d*)', FORMULA)
    atoms = []
    for elem, count in elements:
        count = int(count) if count else 1
        atoms.extend([elem] * count)
    return atoms

@lru_cache(maxsize=100000)
def is_valid_cached(graph_hash):
    """Cached validity check using graph hashes"""
    elements, edges = graph_hash
    g = nx.Graph()
    for i, elem in elements:
        g.add_node(i, element=elem)
    for i, j, order in edges:
        g.add_edge(i, j, order=order)
    return is_valid_uncached(g)

def is_valid_uncached(graph):
    """Original validity check without caching"""
    for node in graph.nodes:
        atom = graph.nodes[node]['element']
        valence = sum(data['order'] for _, _, data in graph.edges(node, data=True))
        if valence > max_valence[atom]:
            return False
    return True

def is_valid(graph):
    """Wrapper function that uses caching"""
    graph_hash = graph_to_hash(graph)
    return is_valid_cached(graph_hash)

def graph_to_hash(g):
    """Create unique hashable representation of graph"""
    elements = tuple(sorted((i, g.nodes[i]['element']) for i in g.nodes))
    edges = tuple(sorted((i, j, g[i][j]['order']) for i, j in g.edges))
    return (elements, edges)

def atoms_graph(atoms):
    """Create base graph from atoms"""
    G = nx.Graph()
    for i, atom in enumerate(atoms):
        G.add_node(i, element=atom)
    return G

def expand_bond_orders(graph):
    """Generate all possible bond order combinations"""
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

def generate_graphs_lazy_worker(args):
    """Worker function for parallel processing"""
    edge_comb, atoms = args
    G_base = atoms_graph(atoms)
    g = G_base.copy()
    valid = True

    for (i, j) in edge_comb:
        a1, a2 = g.nodes[i]['element'], g.nodes[j]['element']
        if (a1, a2) in bond_orders or (a2, a1) in bond_orders:
            g.add_edge(i, j, order=1)
        else:
            valid = False
            break

    results = []
    if valid and nx.is_tree(g) and is_valid(g):
        expanded = expand_bond_orders(g)
        for eg in expanded:
            if is_valid(eg):
                results.append(eg)
    return results

def generate_graphs_lazy(atoms, num_processes=1):
    """Generate molecular graphs with optional multiprocessing"""
    G_base = atoms_graph(atoms)
    n = len(atoms)
    possible_edges = list(combinations(range(n), 2))
    total_combinations = math.comb(len(possible_edges), n - 1)
    
    if num_processes <= 1:
        # Single-process version
        for edge_comb in tqdm(combinations(possible_edges, n - 1), 
                          total=total_combinations, 
                          desc="Generating graphs"):
            for result in generate_graphs_lazy_worker((edge_comb, atoms)):
                yield result
    else:
        # Multiprocessing version
        chunk_size = max(1, total_combinations // (num_processes * 10))
        with Pool(num_processes) as pool:
            for results in tqdm(pool.imap_unordered(
                generate_graphs_lazy_worker,
                ((comb, atoms) for comb in combinations(possible_edges, n - 1)),
                chunksize=chunk_size,
                total=total_combinations
            ), total=total_combinations, desc="Generating graphs"):
                for result in results:
                    yield result

def number_to_bondtype(order):
    """Convert bond order number to RDKit bond type"""
    if order == 1: return Chem.BondType.SINGLE
    if order == 2: return Chem.BondType.DOUBLE
    if order == 3: return Chem.BondType.TRIPLE
    raise ValueError("Invalid bond order")

def graph_to_rdkit_mol(graph):
    """Convert NetworkX graph to RDKit molecule"""
    rw_mol = Chem.RWMol()
    node_to_idx = {}
    for node in graph.nodes:
        atom = Chem.Atom(graph.nodes[node]['element'])
        idx = rw_mol.AddAtom(atom)
        node_to_idx[node] = idx
    for i, j, data in graph.edges(data=True):
        rw_mol.AddBond(node_to_idx[i], node_to_idx[j], number_to_bondtype(data['order']))
    return rw_mol

def generate_smiles(FORMULA, CHARGE, output_file=None, num_processes=1):
    """Main function to generate SMILES with multiprocessing support"""
    output_dir = f"OutputFiles_{FORMULA}_Charge_{CHARGE}"
    os.makedirs(output_dir, exist_ok=True)
    
    if output_file is None:
        output_file = f"nSMILES_{FORMULA}.txt"
    output_path = os.path.join(output_dir, output_file)
    
    atoms = parse_formula(FORMULA)
    graphs = generate_graphs_lazy(atoms, num_processes)
    smiles_set = set()

    for g in tqdm(graphs, desc="Converting to SMILES"):
        mol = graph_to_rdkit_mol(g)
        try:
            smiles = Chem.MolToSmiles(mol, canonical=True)
            smiles_set.add(smiles)
        except:
            continue
    
    with open(output_path, "w") as f:
        for smiles in sorted(smiles_set):
            f.write(smiles + "\n")

    print(f"Saved {len(smiles_set)} SMILES to '{output_path}'")