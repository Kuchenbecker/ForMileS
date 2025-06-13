#########################################################################
# This module is the heart of the program. It generates all graphs      #
# using graph theory and Ullmann's Alogorithm within NetworkX.          #
# Atoms are nodes and bonds are arests. This module generate the list   #
# will be used by all the other modules.                                #
#########################################################################

import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import math
import os

from rdkit import Chem
from rdkit.Chem import Descriptors
from itertools import combinations, product
from rdkit.Chem import Draw
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

def generate_graphs_lazy(atoms):
    G_base = atoms_graph(atoms)
    n = len(atoms)
    possible_edges = list(combinations(range(n), 2))
    total_combinations = math.comb(len(possible_edges), n - 1)

    for edge_comb in tqdm(combinations(possible_edges, n - 1), total=total_combinations, desc="Generating graphs"):
        g = G_base.copy()
        valid = True

        for (i, j) in edge_comb:
            a1, a2 = g.nodes[i]['element'], g.nodes[j]['element']
            if (a1, a2) in bond_orders or (a2, a1) in bond_orders:
                g.add_edge(i, j, order=1)
            else:
                valid = False
                break

        # Avoid cyclic graphs and prune early
        if valid and nx.is_tree(g) and is_valid(g):
            expanded = expand_bond_orders(g)
            for eg in expanded:
                if is_valid(eg):
                    yield eg  


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

def generate_smiles(FORMULA, CHARGE, output_file=None):
    output_dir = f"OutputFiles_{FORMULA}_Charge_{CHARGE}"
    output_path = os.path.join(output_dir, output_file)

    if output_file is None:
        output_file = f"nSMILES_{FORMULA}.txt"
    atoms = parse_formula(FORMULA)
    graphs = generate_graphs_lazy(atoms)
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