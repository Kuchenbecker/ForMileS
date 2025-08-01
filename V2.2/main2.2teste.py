import os
import re
import json
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem import RWMol
from rdkit.Chem.rdchem import BondType
from rdkit.Chem import AllChem
from PIL import Image, ImageDraw, ImageFont
from tqdm import tqdm

# Load configuration
with open("config.json") as f:
    config = json.load(f)

FORMULA = config["FORMULA"]
INPUT_SMILES_LIST = config["INPUT_SMILES_LIST"]
CHARGE = config["CHARGE"]
TARGET_MASS = config["TARGET_MASS"]
TOLERANCE = config["TOLERANCE"]
PARAM_FILE = config["PARAM_FILE"]
SAVE_XYZ = config["SAVE_XYZ"]
SAVE_MOL = config["SAVE_MOL"]
SAVE_SVG = config.get("SAVE_SVG", False)
OUTPUT_DIR = f"OutputFiles_{FORMULA}_Charge_{CHARGE}"

IMG_SIZE = (300, 200)
ANNOTATION_HEIGHT = 60
FONT_SIZE = 14

# Load bond rules
with open(PARAM_FILE, "r") as f:
    params = json.load(f)

bond_orders = {tuple(json.loads(k)): v for k, v in params["bond_orders"].items()}
max_valence = params["max_valence"]
charge_elements = params["charge_elements"]

# Utilities
def create_output_folder():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    if SAVE_XYZ:
        os.makedirs(os.path.join(OUTPUT_DIR, "Coordinate_Files"), exist_ok=True)
    if SAVE_MOL:
        os.makedirs(os.path.join(OUTPUT_DIR, "Mol_Files"), exist_ok=True)

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

# B&B engine
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

# Run generation from multiple SMILES
def run_generation():
    print("[INFO] Starting structure expansion from multiple SMILES...")
    target = parse_formula(FORMULA)
    smiles_list = INPUT_SMILES_LIST

    all_collected = set()

    for smi in tqdm(smiles_list, desc="[GEN] Processing input SMILES"):
        base_mol = Chem.MolFromSmiles(smi)
        if base_mol is None:
            print(f"[WARN] Invalid SMILES: {smi}")
            continue

        rw_base = RWMol(base_mol)
        collected = []
        seen = set()
        grow_recursive(rw_base.GetMol(), target, collected, seen)
        all_collected.update(collected)

    output_path = os.path.join(OUTPUT_DIR, f"nSMILES_{FORMULA}.txt")
    with open(output_path, "w") as f:
        for s in sorted(all_collected):
            f.write(s + "\n")
    print(f"[OK] Generated {len(all_collected)} unique structures.")
    return list(all_collected)

# Charge and mass filtering
def generate_charged_smiles(smiles_list):
    charged = []
    for smi in tqdm(smiles_list, desc="[CHARGE] Adding charges"):
        mol = Chem.MolFromSmiles(smi, sanitize=False)
        if not mol:
            continue
        for atom in mol.GetAtoms():
            if atom.GetSymbol() in charge_elements:
                mol_copy = RWMol(mol)
                atom_idx = atom.GetIdx()
                atom_copy = mol_copy.GetAtomWithIdx(atom_idx)
                atom_copy.SetFormalCharge(CHARGE)
                try:
                    Chem.SanitizeMol(mol_copy)
                    charged_smi = Chem.MolToSmiles(mol_copy, canonical=True)
                    charged.append(charged_smi)
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

# Save molecule images
def smiles_to_images(smiles_list):
    for idx, smiles in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(smiles)
        if not mol: continue

        formula = CalcMolFormula(mol)
        mass = f"{Descriptors.ExactMolWt(mol):.4f}"

        if SAVE_SVG:
            try:
                drawer = Draw.MolDraw2DSVG(IMG_SIZE[0], IMG_SIZE[1])
                drawer.DrawMolecule(mol)
                drawer.FinishDrawing()
                svg = drawer.GetDrawingText()
                fname = f"mol_{idx + 1}.svg"
                with open(os.path.join(OUTPUT_DIR, fname), "w", encoding="utf-8") as f:
                    f.write(svg)
            except Exception as e:
                print(f"[WARN] Failed to generate SVG for molecule {idx+1}: {e}")
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

            fname = f"mol_{idx + 1}.png"
            canvas.save(os.path.join(OUTPUT_DIR, fname))

        if SAVE_XYZ:
            xyz_path = os.path.join(OUTPUT_DIR, "Coordinate_Files", f"mol_{idx + 1}.xyz")
            try:
                AllChem.EmbedMolecule(mol)
                with open(xyz_path, "w") as f:
                    f.write(Chem.MolToXYZBlock(mol))
            except:
                pass

        if SAVE_MOL:
            mol_path = os.path.join(OUTPUT_DIR, "Mol_Files", f"mol_{idx + 1}.mol")
            try:
                with open(mol_path, "w") as f:
                    f.write(Chem.MolToMolBlock(mol))
            except:
                pass

# Main execution
if __name__ == "__main__":
    print("================== ForMileS-SMART v2.2 ==================")
    create_output_folder()
    base_smiles = run_generation()
    charged = generate_charged_smiles(base_smiles)
    final = filter_by_mass(charged)
    smiles_to_images(final)
    print("===================== EXECUTION DONE ====================")
