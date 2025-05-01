import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import math

from rdkit import Chem
from rdkit.Chem import Descriptors
from itertools import combinations, product
from rdkit.Chem import Draw
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from tqdm import tqdm 


# ---- CONFIG ----
CHARGED_TARGET_MASS = 117.092
TARGET_MASS = 118.099
FORMULA = "C6O2"

# ---- Bond and Valence Rules ----
bond_orders = {
    ("C", "C"): [1, 2],
    ("C", "O"): [1, 2],
    ("O", "O"): [1],
}

max_valence = {
    "C": 4,
    "O": 2,
}

# ---- Helper Functions ----

def parse_formula(formula):
    import re
    elements = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
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

def generate_graphs(atoms):
    G_base = atoms_graph(atoms)
    n = len(atoms)
    graphs = []

    # All possible acyclic edge sets with (n - 1) edges
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
                    graphs.append(eg)

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

def generate_smiles(formula):
    atoms = parse_formula(formula)
    graphs = generate_graphs(atoms)
    smiles_list = []

    for g in tqdm(graphs, desc="Converting to SMILES"):
        mol = graph_to_rdkit_mol(g)
        try:
            smiles = Chem.MolToSmiles(mol, canonical=True)
            smiles_list.append(smiles)
        except:
            continue
    return list(set(smiles_list))


def filter_smiles_by_mass(smiles_list, target_mass, tolerance=0.001):
    filtered = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            mass = Descriptors.ExactMolWt(mol)
            if abs(mass - target_mass) <= tolerance:
                filtered.append(smi)
    return filtered

def generate_charged_smiles(smiles_list):
    charged_smiles = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi, sanitize=False)
        if mol is None:
            continue
        for atom in mol.GetAtoms():
            if atom.GetSymbol() in ["C", "O"]:
                mol_copy = Chem.RWMol(mol)
                atom_copy = mol_copy.GetAtomWithIdx(atom.GetIdx())
                atom_copy.SetFormalCharge(1)
                try:
                    charged_smi = Chem.MolToSmiles(mol_copy, canonical=True)
                    charged_smiles.append(charged_smi)
                except:
                    continue
    return list(set(charged_smiles))

def filter_charged_smiles_by_mass(smiles_list, target_mass, tolerance=0.005):
    filtered_charged = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            mass = Descriptors.ExactMolWt(mol)
            if abs(mass - target_mass) <= tolerance:
                filtered_charged.append(smi)
    return filtered_charged

def smiles_to_molecules(smiles_list):
    """Convert SMILES to molecules with mass, formula, and hydrogen count."""
    molecules = []
    data = []  # Will store (smiles, mass, formula)

    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            # Calculate properties
            mass = f"{Descriptors.ExactMolWt(mol):.4f}"
            formula = CalcMolFormula(mol)
            data.append((smiles, mass, formula))
            molecules.append(mol)
        else:
            print(f"Warning: Could not parse SMILES '{smiles}'")

    return molecules, data

def plot_molecules(molecules, data, mols_per_row=4, figsize=(20, 25)):
    """Plot molecules with SMILES, mass, and molecular formula."""
    if not molecules:
        print("No valid molecules to display.")
        return

    n_mols = len(molecules)
    n_rows = math.ceil(n_mols / mols_per_row)

    fig, axes = plt.subplots(n_rows, mols_per_row, figsize=figsize)
    if n_rows == 1:
        axes = axes.reshape(1, -1)

    # Turn off all axes initially
    for ax in axes.flatten():
        ax.axis('off')

    # Plot each molecule
    for i, (mol, (smiles, mass, formula)) in enumerate(zip(molecules, data)):
        row = i // mols_per_row
        col = i % mols_per_row

        # Draw molecule
        img = Draw.MolToImage(mol, size=(300, 200))
        axes[row, col].imshow(img)
        axes[row, col].axis('off')

        # Create annotation text
        annotation = f"{smiles}\nMass: {mass}\n{formula}"

        # Add text annotations
        axes[row, col].text(0.5, 0.01,
                          annotation,
                          transform=axes[row, col].transAxes,
                          ha='center', va='bottom',
                          fontsize=8,
                          bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))

    plt.tight_layout()
    plt.show()


# ---- MAIN PIPELINE ----

if __name__ == "__main__":
    print(f"Generating neutral SMILES for: {FORMULA}")
    smiles_list = generate_smiles(FORMULA)
    print(f"Generated {len(smiles_list)} unique neutral SMILES.\n")

    print(f"Filtering SMILES by exact mass = {TARGET_MASS}")
    filtered = filter_smiles_by_mass(smiles_list, TARGET_MASS)
    print(f"Found {len(filtered)} with mass ≈ {TARGET_MASS}:\n")
    for s in filtered:
        print(s)

    print("\nGenerating charged SMILES from filtered list...")
    charged = generate_charged_smiles(filtered)
    print(f"Generated {len(charged)} charged SMILES:\n")
    for cs in charged:
        print(cs)

    print(f"\nFiltering charged SMILES by exact mass = {CHARGED_TARGET_MASS}")
    filtered_charged = filter_charged_smiles_by_mass(charged, CHARGED_TARGET_MASS)
    print(f"Found {len(filtered_charged)} charged SMILES with mass ≈ {CHARGED_TARGET_MASS}:\n")
    for fcs in filtered_charged:
        print(fcs)

    # Printing the Graphs
    smiles_list = filtered_charged

    # Process molecules and calculate properties
    molecules, data = smiles_to_molecules(smiles_list)

    # Plot the molecules with all information
    plot_molecules(molecules, data, mols_per_row=4)

