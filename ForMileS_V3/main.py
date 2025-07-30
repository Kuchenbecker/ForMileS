#################################################################################
###                                                                           ###
###        ForMileS-SMART v3.1: Multi-SMARTS Structural Expansion Tool        ###
###                                                                           ###
#################################################################################

import os
import re
import json
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw, AllChem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem import RWMol
from rdkit.Chem.rdchem import BondType
from PIL import Image, ImageDraw, ImageFont
from tqdm import tqdm

############################## CONFIGURATION ###################################

FORMULA = "C6O2"
PRECURSOR_SMARTS_LIST = ["CCCOCCCO", "CCCOCCC=O"]
CHARGE = +1
TARGET_MASS = 117.092
TOLERANCE = 0.5
PARAM_FILE = "parameters.json"

SAVE_AS_SVG = True
SAVE_XYZ_FILE = True
SAVE_AS_MOL = True

OUTPUT_DIR = f"OutputFiles_{FORMULA}_Charge_{CHARGE}"
IMG_SIZE = (300, 200)
ANNOTATION_HEIGHT = 60
FONT_SIZE = 14

########################## LOAD BOND/VALENCE RULES #############################
with open(PARAM_FILE, "r") as f:
    params = json.load(f)

bond_orders = {tuple(k): v for k, v in params["bond_orders"].items()}
max_valence = params["max_valence"]
charge_elements = params["charge_elements"]

######################### UTILS & PARSING ######################################
def create_output_folder():
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)

def parse_formula(formula):
    matches = re.findall(r"([A-Z][a-z]*)(\d*)", formula)
    atom_counts = {}
    for elem, count in matches:
        count = int(count) if count else 1
        atom_counts[elem] = atom_counts.get(elem, 0) + count
    return atom_counts

def count_atoms(mol):
    atom_counts = {}
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        atom_counts[symbol] = atom_counts.get(symbol, 0) + 1
    return atom_counts

def atom_deficit(target, current):
    return {k: target.get(k, 0) - current.get(k, 0) for k in target}

def is_deficit_valid(deficit):
    return all(v >= 0 for v in deficit.values())

def valence_ok(mol):
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        valence = sum([b.GetBondTypeAsDouble() for b in atom.GetBonds()])
        if valence > max_valence.get(symbol, 4):
            return False
    return True

def open_sites(mol):
    result = []
    for atom in mol.GetAtoms():
        symbol = atom.GetSymbol()
        if symbol not in max_valence:
            continue
        valence = sum([b.GetBondTypeAsDouble() for b in atom.GetBonds()])
        if valence < max_valence[symbol]:
            result.append(atom.GetIdx())
    return result

def attach_atom(base_mol, atom_symbol, attach_to_idx, bond_order):
    mol = RWMol(base_mol)
    new_atom = Chem.Atom(atom_symbol)
    new_idx = mol.AddAtom(new_atom)
    mol.AddBond(attach_to_idx, new_idx, BondType.values[bond_order])
    return mol.GetMol()

def mol_to_canonical_smiles(mol):
    try:
        return Chem.MolToSmiles(mol, canonical=True)
    except:
        return None

########################### B&B EXPANSION ENGINE ###############################
def grow_recursive(current_mol, target_formula, collected, seen):
    current_atoms = count_atoms(current_mol)
    deficit = atom_deficit(target_formula, current_atoms)

    if not is_deficit_valid(deficit):
        return

    if sum(deficit.values()) == 0:
        smiles = mol_to_canonical_smiles(current_mol)
        if smiles and smiles not in seen:
            seen.add(smiles)
            collected.append(smiles)
        return

    if not valence_ok(current_mol):
        return

    for atom_symbol in deficit:
        if deficit[atom_symbol] == 0:
            continue
        for idx in open_sites(current_mol):
            symbol1 = current_mol.GetAtomWithIdx(idx).GetSymbol()
            pair = tuple(sorted((symbol1, atom_symbol)))
            for order in bond_orders.get(pair, []):
                try:
                    new_mol = attach_atom(current_mol, atom_symbol, idx, order)
                    grow_recursive(new_mol, target_formula, collected, seen)
                except:
                    continue

############################## MOLECULE GENERATION #############################
def run_generation():
    print("[INFO] Starting SMARTS expansion...")
    target = parse_formula(FORMULA)
    collected = []
    seen = set()

    for smarts in PRECURSOR_SMARTS_LIST:
        base_mol = Chem.MolFromSmarts(smarts)
        if base_mol is None:
            print(f"[WARNING] Could not parse SMARTS: {smarts}")
            continue
        rw_base = RWMol(base_mol)
        base_mol = rw_base.GetMol()
        grow_recursive(base_mol, target, collected, seen)

    output_path = os.path.join(OUTPUT_DIR, f"nSMILES_{FORMULA}.txt")
    with open(output_path, "w") as f:
        for s in sorted(collected):
            f.write(s + "\n")

    print(f"[OK] Saved {len(collected)} SMILES to {output_path}")
    return collected

############################## ADD CHARGE TO ATOMS #############################
def generate_charged_smiles(smiles_list):
    charged = []
    for smi in tqdm(smiles_list, desc="[CHARGE] Adding charge"):
        mol = Chem.MolFromSmiles(smi)
        if not mol: continue
        for atom in mol.GetAtoms():
            if atom.GetSymbol() in charge_elements:
                mol_copy = RWMol(mol)
                mol_copy.GetAtomWithIdx(atom.GetIdx()).SetFormalCharge(CHARGE)
                try:
                    csmi = Chem.MolToSmiles(mol_copy, canonical=True)
                    charged.append(csmi)
                except:
                    continue
    path = os.path.join(OUTPUT_DIR, f"chargedSMILES_{FORMULA}.txt")
    with open(path, "w") as f:
        for c in charged:
            f.write(c + "\n")
    return charged

def filter_by_mass(smiles_list):
    filtered = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if not mol: continue
        mass = Descriptors.ExactMolWt(mol)
        if abs(mass - TARGET_MASS) <= TOLERANCE:
            filtered.append(smi)
    path = os.path.join(OUTPUT_DIR, f"filteredchargedSMILES_{FORMULA}.txt")
    with open(path, "w") as f:
        for c in filtered:
            f.write(c + "\n")
    return filtered

############################## VISUALIZATION + XYZ #############################
def sanitize_filename(smiles):
    return re.sub(r'[^a-zA-Z0-9._-]', '_', smiles)

def smiles_to_images_and_xyz(smiles_list):
    coord_dir = os.path.join(OUTPUT_DIR, "coordinate_files")
    mol_dir = os.path.join(OUTPUT_DIR, "mol_files")

    if SAVE_XYZ_FILE:
        os.makedirs(coord_dir, exist_ok=True)
    if SAVE_AS_MOL:
        os.makedirs(mol_dir, exist_ok=True)

    for idx, smiles in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(smiles)
        if not mol: continue

        formula = CalcMolFormula(mol)
        mass = f"{Descriptors.ExactMolWt(mol):.4f}"
        fname_base = f"mol_{idx + 1}"

        # Image output
        if SAVE_AS_SVG:
            svg_path = os.path.join(OUTPUT_DIR, fname_base + ".svg")
            try:
                Draw.MolToFile(mol, svg_path, size=IMG_SIZE, imageType="svg")
            except Exception as e:
                print(f"[WARNING] Failed SVG render for {smiles}: {e}")
        else:
            img = Draw.MolToImage(mol, size=IMG_SIZE)
            total_height = IMG_SIZE[1] + ANNOTATION_HEIGHT
            canvas = Image.new("RGB", (IMG_SIZE[0], total_height), "white")
            canvas.paste(img, (0, 0))
            draw = ImageDraw.Draw(canvas)
            try:
                font = ImageFont.truetype("arial.ttf", FONT_SIZE)
            except:
                font = ImageFont.load_default()
            draw.text((10, IMG_SIZE[1] + 5), f"{formula} | {mass}", fill="black", font=font)
            draw.text((10, IMG_SIZE[1] + 25), smiles, fill="black", font=font)
            img_path = os.path.join(OUTPUT_DIR, fname_base + ".png")
            canvas.save(img_path)

        # XYZ and MOL coordinate files
        if SAVE_XYZ_FILE or SAVE_AS_MOL:
            mol_with_H = Chem.AddHs(mol)
            try:
                AllChem.EmbedMolecule(mol_with_H, AllChem.ETKDG())
                conf = mol_with_H.GetConformer()

                if SAVE_XYZ_FILE:
                    xyz_path = os.path.join(coord_dir, fname_base + ".xyz")
                    with open(xyz_path, "w") as f:
                        f.write(f"{mol_with_H.GetNumAtoms()}\n{smiles}\n")
                        for atom in mol_with_H.GetAtoms():
                            pos = conf.GetAtomPosition(atom.GetIdx())
                            f.write(f"{atom.GetSymbol():<2} {pos.x:.4f} {pos.y:.4f} {pos.z:.4f}\n")

                if SAVE_AS_MOL:
                    mol_path = os.path.join(mol_dir, fname_base + ".mol")
                    Chem.MolToMolFile(mol_with_H, mol_path)

            except:
                print(f"[WARNING] Could not generate 3D coordinates for: {smiles}")

############################## EXECUTION #######################################
if __name__ == "__main__":
    print("================== ForMileS v3.1 ==================")
    create_output_folder()
    base_smiles = run_generation()
    charged = generate_charged_smiles(base_smiles)
    final = filter_by_mass(charged)
    smiles_to_images_and_xyz(final)
    print("=================== EXECUTION DONE ===================")
