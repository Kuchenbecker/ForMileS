
import os
import re
import json
import math
import networkx as nx
from itertools import combinations, product
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem import RWMol
from rdkit.Chem.rdchem import BondType
from PIL import Image, ImageDraw, ImageFont
from tqdm import tqdm

# CONFIGURAÇÕES
FORMULA = "C6O2"
CHARGE = 1
TARGET_MASS = 117.103
TOLERANCE = 0.5
OUTPUT_DIR = "OutputFiles_" + FORMULA + "_Charge_" + str(CHARGE)
IMG_SIZE = (300, 200)
ANNOTATION_HEIGHT = 60
FONT_SIZE = 14
SAVE_XYZ = False
SAVE_MOL = False
PRECURSOR_SMARTS = "CCCOCCCO"

bond_orders = {
    ('C', 'C'): [1, 2, 3],
    ('C', 'O'): [1, 2],
    ('O', 'O'): [1],
    ('C', 'N'): [1, 2, 3],
    ('N', 'O'): [1]
}
max_valence = {'C': 4, 'O': 2, 'N': 3}
charge_elements = ['O', 'N', 'C']

def create_output_folder():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

def parse_formula(formula):
    matches = re.findall(r"([A-Z][a-z]*)(\d*)", formula)
    atom_list = []
    for elem, count in matches:
        atom_list.extend([elem] * (int(count) if count else 1))
    return atom_list

def atoms_graph(atoms):
    g = nx.Graph()
    for i, atom in enumerate(atoms):
        g.add_node(i, element=atom)
    return g

def is_valid(g):
    valences = dict((n, 0) for n in g.nodes)
    for i, j, data in g.edges(data=True):
        valences[i] += data.get("order", 1)
        valences[j] += data.get("order", 1)
    for i in g.nodes:
        symbol = g.nodes[i]["element"]
        if valences[i] > max_valence.get(symbol, 4):
            return False
    return True

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

def generate_neutral_smiles():
    atoms = parse_formula(FORMULA)
    G_base = atoms_graph(atoms)
    n = len(atoms)
    possible_edges = list(combinations(range(n), 2))
    smiles_set = set()
    for edge_comb in combinations(possible_edges, n - 1):
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
                    mol = graph_to_rdkit_mol(eg)
                    try:
                        smi = Chem.MolToSmiles(mol, canonical=True)
                        smiles_set.add(smi)
                    except:
                        continue
    path = os.path.join(OUTPUT_DIR, f"nSMILES_" + FORMULA + ".txt")
    with open(path, "w") as f:
        for smi in sorted(smiles_set):
            f.write(smi + "\n")
    return list(smiles_set)

def generate_charged_smiles(smiles_list):
    charged = []
    for smi in tqdm(smiles_list, desc=f"[CHARGE] Adicionando carga"):
        mol = Chem.MolFromSmiles(smi, sanitize=False)
        if not mol:
            continue
        charge_sites = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetSymbol() in charge_elements]
        if len(charge_sites) < CHARGE:
            continue
        for combo in combinations(charge_sites, CHARGE):
            mol_copy = RWMol(mol)
            for idx in combo:
                mol_copy.GetAtomWithIdx(idx).SetFormalCharge(1)
            try:
                charged_smi = Chem.MolToSmiles(mol_copy, canonical=True)
                charged.append(charged_smi)
            except:
                continue
    path = os.path.join(OUTPUT_DIR, f"chargedSMILES_" + FORMULA + ".txt")
    with open(path, "w") as f:
        for smi in sorted(set(charged)):
            f.write(smi + "\n")
    return charged

def filter_by_mass(smiles_list):
    filtered = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if not mol: continue
        mass = Descriptors.ExactMolWt(mol)
        if abs(mass - TARGET_MASS) <= TOLERANCE:
            filtered.append(smi)
    path = os.path.join(OUTPUT_DIR, f"filteredchargedSMILES_" + FORMULA + ".txt")
    with open(path, "w") as f:
        for smi in filtered:
            f.write(smi + "\n")
    return filtered

def smiles_to_images(smiles_list):
    for idx, smiles in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(smiles)
        if not mol: continue
        formula = CalcMolFormula(mol)
        img = Draw.MolToImage(mol, size=IMG_SIZE)
        annotated_img = Image.new("RGB", (IMG_SIZE[0], IMG_SIZE[1] + ANNOTATION_HEIGHT), "white")
        annotated_img.paste(img, (0, 0))
        draw = ImageDraw.Draw(annotated_img)
        font = ImageFont.load_default()
        draw.text((5, IMG_SIZE[1] + 5), formula, fill="black", font=font)
        img_path = os.path.join(OUTPUT_DIR, f"mol_" + str(idx + 1) + ".png")
        annotated_img.save(img_path)

if __name__ == "__main__":
    create_output_folder()
    neutral = generate_neutral_smiles()
    charged = generate_charged_smiles(neutral)
    filtered = filter_by_mass(charged)
    smiles_to_images(filtered)
