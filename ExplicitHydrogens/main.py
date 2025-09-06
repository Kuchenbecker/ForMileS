#################################################################################
###                                                                           ###
###                   ForMileS: Formation of Mass SMILES                      ###
###                     v2.6 — Explicit Hydrogens Ready                       ###
###                                                                           ###
#################################################################################

import os
import re
import json
import time
import tracemalloc
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from rdkit.Chem import RWMol
from rdkit.Chem.rdchem import BondType
from rdkit.Chem import AllChem
from PIL import Image, ImageDraw, ImageFont
from tqdm import tqdm
import sys

def resource_path(relative_path):
    """ Get absolute path to resource, works for dev and for PyInstaller """
    try:
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.abspath(".")

    return os.path.join(base_path, relative_path)

############################## LOAD CONFIG ######################################
with open(resource_path("config.json")) as f:
    config = json.load(f)

FORMULA = config["FORMULA"]
MOL_SCAFFOLD = config["MOL_SCAFFOLD"]
CHARGE = config["CHARGE"]
TARGET_MASS = config["TARGET_MASS"]
TOLERANCE = config["TOLERANCE"]
PARAM_FILE = config["PARAM_FILE"]
SAVE_XYZ = config["SAVE_XYZ"]
SAVE_MOL = config["SAVE_MOL"]
ALLOW_BRANCHING = config.get("ALLOW_BRANCHING")
ALLOW_CYCLES = config.get("ALLOW_CYCLES")
ALLOW_DOUBLE_BONDS = config.get("ALLOW_DOUBLE_BONDS", True)
ALLOW_TRIPLE_BONDS = config.get("ALLOW_TRIPLE_BONDS", False)
MAX_DOUBLE_BONDS = config.get("MAX_DOUBLE_BONDS", 2)
MAX_TRIPLE_BONDS = config.get("MAX_TRIPLE_BONDS", 0)  # Default to no triple bonds
OUTPUT_DIR = config.get("OUTPUT_DIR", f"OutputFiles_{FORMULA}_Charge_{CHARGE}")
EXPLICIT_HYDROGENS = config.get("EXPLICIT_HYDROGENS", False)
BENCHMARK = config.get("BENCHMARK", False)

IMG_SIZE = (300, 200)
ANNOTATION_HEIGHT = 60
FONT_SIZE = 14

############## LOAD BOND AND VALENCE RULES FROM FILE ###########################
with open(resource_path(config["PARAM_FILE"]), "r") as f:  # Modified line
    params = json.load(f)

bond_orders = {tuple(json.loads(k)): v for k, v in params["bond_orders"].items()}
max_valence = params["max_valence"]
charge_elements = params["charge_elements"]

######################### UTILITY FUNCTIONS ####################################
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
        valence = atom.GetTotalValence()  # More reliable method
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
        # Preserve explicit Hs if the mode is ON
        return Chem.MolToSmiles(mol, canonical=True, allHsExplicit=EXPLICIT_HYDROGENS)
    except:
        return None

def is_linear(mol):
    """Check if molecule is strictly linear (no branching)"""
    for atom in mol.GetAtoms():
        if len(atom.GetBonds()) > 2:
            return False
    return True

def has_cycles(mol):
    """Check if molecule contains any rings"""
    return mol.GetRingInfo().NumRings() > 0

def count_bond_types(mol):
    """Count number of double and triple bonds in molecule"""
    dbl = 0
    trpl = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == BondType.DOUBLE:
            dbl += 1
        elif bond.GetBondType() == BondType.TRIPLE:
            trpl += 1
    return dbl, trpl

def get_allowed_bond_orders(pair):
    """Filter bond orders based on configuration"""
    orders = bond_orders.get(pair, [1])  # Default to single bond if pair not defined
    filtered = []
    for order in orders:
        if order == 1:
            filtered.append(order)  # Always allow single bonds
        elif order == 2 and ALLOW_DOUBLE_BONDS and MAX_DOUBLE_BONDS > 0:
            filtered.append(order)
        elif order == 3 and ALLOW_TRIPLE_BONDS and MAX_TRIPLE_BONDS > 0:
            filtered.append(order)
    return filtered

####################### B&B EXPANSION ENGINE ####################################
def explore_bond_permutations(mol, target_formula, collected, seen, depth=0, max_depth=3):
    """Explore different bond order combinations for completed molecules"""
    if depth >= max_depth:  # Prevent infinite recursion
        return

    rw_mol = RWMol(mol)
    bonds = list(rw_mol.GetBonds())
    current_dbl, current_trpl = count_bond_types(rw_mol)

    for bond in bonds:
        a1 = bond.GetBeginAtom().GetSymbol()
        a2 = bond.GetEndAtom().GetSymbol()
        pair = tuple(sorted((a1, a2)))
        current_order = bond.GetBondTypeAsDouble()

        # Skip if current bond order is already maximum
        if current_order == 3:
            continue

        # Get allowed higher bond orders for this pair
        possible_orders = [o for o in get_allowed_bond_orders(pair) if o > current_order]

        for new_order in possible_orders:
            # Skip if we'd exceed max allowed bonds
            if new_order == 2 and current_dbl >= MAX_DOUBLE_BONDS:
                continue
            if new_order == 3 and current_trpl >= MAX_TRIPLE_BONDS:
                continue

            # Create new molecule with modified bond
            new_mol = RWMol(rw_mol)
            new_bond = new_mol.GetBondWithIdx(bond.GetIdx())
            new_bond.SetBondType(BondType.values[new_order])

            try:
                Chem.SanitizeMol(new_mol)
                process_complete_molecule(new_mol, target_formula, collected, seen, depth+1)
            except:
                continue

def process_complete_molecule(mol, target_formula, collected, seen, depth=0):
    """Handle molecules that match the target formula"""
    smiles = mol_to_canonical_smiles(mol)
    if not smiles or smiles in seen:
        return

    # Check structural constraints
    if (not ALLOW_BRANCHING and not is_linear(mol)) or (not ALLOW_CYCLES and has_cycles(mol)):
        return

    seen.add(smiles)
    collected.append(smiles)

    # Explore bond permutations if allowed
    if (ALLOW_DOUBLE_BONDS and MAX_DOUBLE_BONDS > 0) or (ALLOW_TRIPLE_BONDS and MAX_TRIPLE_BONDS > 0):
        explore_bond_permutations(mol, target_formula, collected, seen, depth)

def _atoms_with_replaceable_h(mol):
    """Return list of (atom_idx, h_neighbor_idx) where atom has at least one attached hydrogen.
    Only used when EXPLICIT_HYDROGENS=True to create growth sites by swapping H for a heavy atom.
    """
    pairs = []
    for a in mol.GetAtoms():
        for nb in a.GetNeighbors():
            if nb.GetSymbol() == "H":
                pairs.append((a.GetIdx(), nb.GetIdx()))
                break  # one H is enough
    return pairs

def _swap_h_for_atom(base_mol, host_idx, h_idx, new_symbol, bond_order=1):
    """Remove hydrogen (h_idx) attached to host(atom), and add a new atom (new_symbol) with the same bond position."""
    rw = RWMol(base_mol)
    # Remove the H first (ensure removing the correct atom index by ordering)
    # Always delete the higher index first to avoid reindexing issues
    del_h_idx = h_idx
    rw.RemoveAtom(del_h_idx)
    # After removal, host index may shift if h_idx < host_idx
    if del_h_idx < host_idx:
        host_idx -= 1
    # Add new atom and bond
    new_idx = rw.AddAtom(Chem.Atom(new_symbol))
    rw.AddBond(host_idx, new_idx, BondType.values[bond_order])
    return rw.GetMol()

def grow_recursive(current_mol, target_formula, collected, seen, depth=0):
    current_atoms = count_atoms(current_mol)
    deficit = atom_deficit(target_formula, current_atoms)

    # BOUNDING STEP 1: Prune if atom deficit is invalid
    if not is_deficit_valid(deficit):
        return

    # BOUNDING STEP 2: Check if molecule is complete
    if sum(deficit.values()) == 0:
        process_complete_molecule(current_mol, target_formula, collected, seen, depth)
        return

    # BOUNDING STEP 3: Prune if valence rules are violated
    if not valence_ok(current_mol):
        return

    # BOUNDING STEP 4: Prune based on structural constraints
    if (not ALLOW_BRANCHING and not is_linear(current_mol)) or (not ALLOW_CYCLES and has_cycles(current_mol)):
        return

    # Determine open growth sites
    sites = open_sites(current_mol)

    # If there are no open sites AND explicit-H mode is on, try "H-swap" to create a site:
    # replace one attached hydrogen on a saturated atom by the new heavy atom to grow the graph.
    if EXPLICIT_HYDROGENS and (len(sites) == 0):
        # Consider only heavy atoms (exclude H) for swapping-in
        heavy_symbols = [s for s in deficit if s != "H" and deficit[s] > 0]
        if heavy_symbols:
            repl_candidates = _atoms_with_replaceable_h(current_mol)
            for host_idx, h_idx in repl_candidates:
                symbol1 = current_mol.GetAtomWithIdx(host_idx).GetSymbol()
                # Try each heavy element still missing
                for atom_symbol in heavy_symbols:
                    pair = tuple(sorted((symbol1, atom_symbol)))
                    for order in get_allowed_bond_orders(pair):
                        try:
                            new_mol = _swap_h_for_atom(current_mol, host_idx, h_idx, atom_symbol, bond_order=order)
                            Chem.SanitizeMol(new_mol)
                            grow_recursive(new_mol, target_formula, collected, seen, depth+1)
                        except Exception:
                            continue
        return  # either we branched via swaps or there's nothing to do

    # BRANCHING: Explore all possible atom additions (including H when EXPLICIT_HYDROGENS and present in target)
    for atom_symbol in deficit:
        if deficit[atom_symbol] == 0:
            continue
        for idx in sites:
            symbol1 = current_mol.GetAtomWithIdx(idx).GetSymbol()
            pair = tuple(sorted((symbol1, atom_symbol)))
            for order in get_allowed_bond_orders(pair):
                try:
                    new_mol = attach_atom(current_mol, atom_symbol, idx, order)
                    Chem.SanitizeMol(new_mol)  # Sanitize after adding atom
                    grow_recursive(new_mol, target_formula, collected, seen, depth+1)
                except Exception:
                    continue

    current_atoms = count_atoms(current_mol)
    deficit = atom_deficit(target_formula, current_atoms)

    # BOUNDING STEP 1: Prune if atom deficit is invalid
    if not is_deficit_valid(deficit):
        return

    # BOUNDING STEP 2: Check if molecule is complete
    if sum(deficit.values()) == 0:
        process_complete_molecule(current_mol, target_formula, collected, seen, depth)
        return

    # BOUNDING STEP 3: Prune if valence rules are violated
    if not valence_ok(current_mol):
        return

    # BOUNDING STEP 4: Prune based on structural constraints
    if (not ALLOW_BRANCHING and not is_linear(current_mol)) or (not ALLOW_CYCLES and has_cycles(current_mol)):
        return

    # BRANCHING: Explore all possible atom additions (including H when EXPLICIT_HYDROGENS and present in target)
    for atom_symbol in deficit:
        if deficit[atom_symbol] == 0:
            continue
        for idx in open_sites(current_mol):
            symbol1 = current_mol.GetAtomWithIdx(idx).GetSymbol()
            pair = tuple(sorted((symbol1, atom_symbol)))
            for order in get_allowed_bond_orders(pair):
                try:
                    new_mol = attach_atom(current_mol, atom_symbol, idx, order)
                    Chem.SanitizeMol(new_mol)  # Sanitize after adding atom
                    grow_recursive(new_mol, target_formula, collected, seen, depth+1)
                except:
                    continue

####################### MOL GRAPH GENERATION ####################################
def run_generation():
    if BENCHMARK:
        start_time = time.perf_counter()
        tracemalloc.start()

    print("[INFO] Initiating molecular graph expansion from scaffolds...")
    print(f"[MODE] Explicit hydrogens: {EXPLICIT_HYDROGENS}")
    target = parse_formula(FORMULA)

    if isinstance(MOL_SCAFFOLD, str):
        precursor_list = [MOL_SCAFFOLD]
    else:
        precursor_list = MOL_SCAFFOLD

    collected = []
    seen = set()

    for smi in precursor_list:
        print(f"[PROCESSING] Building from scaffold: {smi}")
        try:
            # Preserve explicit Hs on input SMILES if present
            base_mol = Chem.MolFromSmiles(smi, sanitize=True)
            if not base_mol:
                print(f"[WARNING] Invalid SMILES: {smi}")
                continue
            if EXPLICIT_HYDROGENS:
                # Materialize ALL implicit hydrogens on the starting scaffold to avoid valence overflows
                # and to enable H-swap growth on saturated atoms.
                base_mol = Chem.AddHs(base_mol, addCoords=False)
            else:
                # Keep implicit-H behavior (no explicit Hs on scaffold)
                # Nothing to do here.
                pass

            # Sanitize the initial molecule
            Chem.SanitizeMol(base_mol)

            rw_base = RWMol(base_mol)
            base_mol = rw_base.GetMol()
            grow_recursive(base_mol, target, collected, seen)
        except Exception as e:
            print(f"[ERROR] Processing scaffold {smi}: {str(e)}")
            continue

    output_path = os.path.join(OUTPUT_DIR, f"nSMILES_{FORMULA}.txt")
    with open(output_path, "w") as f:
        for s in sorted(collected):
            f.write(s + "\n")

    print(f"[OK] Saved {len(collected)} unique SMILES in {output_path}")

    if BENCHMARK:
        runtime_s = time.perf_counter() - start_time
        current, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        with open(os.path.join(OUTPUT_DIR, "benchmark.txt"), "w") as bf:
            bf.write(f"Explicit Hs: {EXPLICIT_HYDROGENS}\n")
            bf.write(f"Runtime (s): {runtime_s:.3f}\n")
            bf.write(f"Peak memory (MiB): {peak/1024/1024:.2f}\n")
        print(f"[BENCH] Runtime={runtime_s:.3f}s | Peak={peak/1024/1024:.2f} MiB")

    return collected

####################### CHARGE AND MASS FILTERING ####################################
def generate_charged_smiles(smiles_list):
    charged = []
    for smi in tqdm(smiles_list, desc="[CHARGE] Adding charge"):
        # Preserve explicit Hs if already present in the SMILES string
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
                    charged_smi = Chem.MolToSmiles(
                        mol_copy, canonical=True, allHsExplicit=EXPLICIT_HYDROGENS
                    )
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

######################### OUTPUT REPORTING ##########################################
def _maybe_add_hs(mol):
    """Return a molecule with hydrogens handled according to mode.
       - If EXPLICIT_HYDROGENS is True: assume Hs explicitly present are the source of truth (do NOT AddHs).
       - Else: expand to explicit Hs for 3D generation/depictions.
    """
    if EXPLICIT_HYDROGENS:
        return Chem.Mol(mol)  # shallow copy, keep as-is
    else:
        return Chem.AddHs(mol)

def generate_xyz_from_smiles(smiles):
    """Generate optimized XYZ coordinates from SMILES string"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"[WARNING] Could not parse SMILES: {smiles}")
            return None

        mol_with_h = _maybe_add_hs(mol)

        # Generate 3D coordinates
        if AllChem.EmbedMolecule(mol_with_h) != 0:
            print(f"[WARNING] 3D embedding failed for: {smiles}")
            return None

        # Optimize geometry (UFF is faster than MMFF for large molecules)
        if AllChem.UFFOptimizeMolecule(mol_with_h) != 0:
            print(f"[WARNING] Optimization failed for: {smiles}")
            return None

        return Chem.MolToXYZBlock(mol_with_h)
    except Exception as e:
        print(f"[ERROR] XYZ generation failed for {smiles}: {str(e)}")
        return None

def generate_mol_from_smiles(smiles):
    """Generate MOL file content from SMILES string"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"[WARNING] Could not parse SMILES: {smiles}")
            return None

        mol_with_h = _maybe_add_hs(mol)

        if AllChem.EmbedMolecule(mol_with_h) != 0:
            print(f"[WARNING] 3D embedding failed for: {smiles}")
            return None

        # Optimize geometry
        if AllChem.UFFOptimizeMolecule(mol_with_h) != 0:
            print(f"[WARNING] Optimization failed for: {smiles}")
            return None

        return Chem.MolToMolBlock(mol_with_h)       
    except Exception as e:
        print(f"[ERROR] MOL generation failed for {smiles}: {str(e)}")
        return None

def smiles_to_images(smiles_list):
    # Create SVG folder if needed
    if config.get("SAVE_SVG", False):
        svg_dir = os.path.join(OUTPUT_DIR, "SVG_Files")
        os.makedirs(svg_dir, exist_ok=True)

    for idx, smiles in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(smiles)
        if not mol: 
            continue

        # Generate 2D coordinates
        try:
            AllChem.Compute2DCoords(mol)
        except:
            print(f"[WARNING] Failed to generate 2D coordinates for {smiles}")
            continue

        # For annotation, respect explicit-H mode in the formula display
        try:
            if EXPLICIT_HYDROGENS:
                # CalcMolFormula uses explicit atoms present;
                # to ensure formula reflects explicit Hs, keep as-is.
                formula = CalcMolFormula(mol)
            else:
                # When not explicit, RDKit formula already includes implicit Hs
                formula = CalcMolFormula(mol)
        except Exception:
            formula = "?"

        mass = f"{Descriptors.ExactMolWt(mol):.4f}"

        # Always generate PNG with annotations
        try:
            # When drawing, show explicit Hs in atom labels if the user wants to see them
            img = Draw.MolToImage(mol, size=IMG_SIZE, kekulize=True, options=None)
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
        except Exception as e:
            print(f"[WARNING] Failed to generate PNG for {smiles}: {str(e)}")
            continue

        # Optionally generate SVG (just the structure)
        if config.get("SAVE_SVG", False):
            try:
                drawer = Draw.MolDraw2DSVG(IMG_SIZE[0], IMG_SIZE[1])
                drawer.DrawMolecule(mol)
                drawer.FinishDrawing()
                svg_data = drawer.GetDrawingText()

                svg_path = os.path.join(OUTPUT_DIR, "SVG_Files", f"mol_{idx + 1}.svg")
                with open(svg_path, "w") as f:
                    f.write(svg_data)
            except Exception as e:
                print(f"[WARNING] Failed to generate SVG for {smiles}: {str(e)}")
                continue

        # Handle XYZ and MOL files
        if SAVE_XYZ:
            xyz_path = os.path.join(OUTPUT_DIR, "Coordinate_Files", f"mol_{idx + 1}.xyz")
            try:
                xyz_block = generate_xyz_from_smiles(smiles)
                if xyz_block:
                    with open(xyz_path, "w") as f:
                        f.write(xyz_block)
            except Exception as e:
                print(f"[WARNING] Failed to generate XYZ for {smiles}: {str(e)}")

        if SAVE_MOL:
            mol_path = os.path.join(OUTPUT_DIR, "Mol_Files", f"mol_{idx + 1}.mol")
            try:
                mol_block = generate_mol_from_smiles(smiles)
                if mol_block:
                    with open(mol_path, "w") as f:
                        f.write(mol_block)
            except Exception as e:
                print(f"[WARNING] Failed to generate MOL for {smiles}: {str(e)}")

########################### MAIN EXECUTION ######################################
if __name__ == "__main__":
    print("================== ForMileS v2.6 (Explicit-H) ==================")
    print(f"Structural Setup: Molecular Branching={ALLOW_BRANCHING}, Cyclic={ALLOW_CYCLES}")
    create_output_folder()
    base_smiles = run_generation()
    charged = generate_charged_smiles(base_smiles)
    final = filter_by_mass(charged)
    smiles_to_images(final)
    print("=================== FINISHED ;) ===================")