#################################################################################
###                                                                           ###
###                   ForMileS: Formation of Mass SMILES                      ###
###                                                                           ###
#################################################################################

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
import sys
import time

# --- Minimal RAM helper (prefers psutil; falls back to Unix resource) ---
def _get_ram_mb():
    try:
        import psutil
        return psutil.Process().memory_info().rss / (1024 * 1024)  # current RSS
    except Exception:
        try:
            import resource, platform
            r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
            # macOS returns bytes; Linux returns kilobytes
            return (r / (1024 * 1024)) if platform.system() == "Darwin" else (r / 1024.0)
        except Exception:
            return None


def resource_path(relative_path):
    """ Get absolute path to resource, works for dev and for PyInstaller """
    try:
        # PyInstaller creates a temp folder and stores path in _MEIPASS
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
        return Chem.MolToSmiles(mol, canonical=True)
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

# ---------------------------- RING CLOSURE (NEW) ------------------------------
def _current_valence(atom):
    """Sum of bond orders around atom."""
    return sum(b.GetBondTypeAsDouble() for b in atom.GetBonds())

def _has_bond(mol, i, j):
    return mol.GetBondBetweenAtoms(i, j) is not None

def try_ring_closures(mol, target_formula, collected, seen, depth=0, max_new=10):

    if not ALLOW_CYCLES:
        return

    new_count = 0
    rw = RWMol(mol)
    n = rw.GetNumAtoms()

    for i in range(n):
        ai = rw.GetAtomWithIdx(i)
        si = ai.GetSymbol()
        if _current_valence(ai) >= max_valence.get(si, 4):
            continue

        for j in range(i + 1, n):
            aj = rw.GetAtomWithIdx(j)
            sj = aj.GetSymbol()

            if _has_bond(rw, i, j):
                continue

            if _current_valence(aj) >= max_valence.get(sj, 4):
                continue

            pair = tuple(sorted((si, sj)))
            possible_orders = bond_orders.get(pair, [1])
            if 1 not in possible_orders:
                continue

            new_mol = RWMol(rw)
            new_mol.AddBond(i, j, BondType.SINGLE)

            try:
                Chem.SanitizeMol(new_mol)
            except Exception:
                continue

            if new_mol.GetRingInfo().NumRings() == 0:
                continue

            # Process as a completed molecule (handles seen/dedup + permutations)
            process_complete_molecule(new_mol, target_formula, collected, seen, depth + 1)

            new_count += 1
            if new_count >= max_new:
                return  # limit the number of closures per parent

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

    # NEW: if cycles are allowed, attempt ring closures from this parent
    if ALLOW_CYCLES:
        try_ring_closures(mol, target_formula, collected, seen, depth)

    # Explore bond permutations if allowed
    if (ALLOW_DOUBLE_BONDS and MAX_DOUBLE_BONDS > 0) or (ALLOW_TRIPLE_BONDS and MAX_TRIPLE_BONDS > 0):
        explore_bond_permutations(mol, target_formula, collected, seen, depth)

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

    # BRANCHING: Explore all possible atom additions
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
    print("[INFO] Initiating molecular graph expansion from scaffolds...")
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
            base_mol = Chem.MolFromSmiles(smi)
            if not base_mol:
                print(f"[WARNING] Invalid SMILES: {smi}")
                continue
                
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
    return collected

####################### CHARGE AND MASS FILTERING ####################################
def generate_charged_smiles(smiles_list):
    """
    Charge exactly one representative atom per symmetry class among chargeable elements,
    and deduplicate by canonical SMILES before returning and writing to disk.
    """
    charged = []
    seen = set()

    for smi in tqdm(smiles_list, desc="[CHARGE] Adding charge"):
        mol = Chem.MolFromSmiles(smi, sanitize=False)
        if not mol:
            continue

        # Symmetry classes for atoms (stable across RDKit versions)
        try:
            ranks = list(Chem.CanonicalRankAtoms(mol))
        except Exception:
            # Fallback: unique rank per atom (disables symmetry pruning if RDKit lacks the call)
            ranks = list(range(mol.GetNumAtoms()))

        # Choose one representative atom index per symmetry rank (only for allowed elements)
        reps = {}
        for a in mol.GetAtoms():
            if a.GetSymbol() in charge_elements:
                r = ranks[a.GetIdx()]
                # Use (element, rank) so different elements with same rank don’t collide
                reps.setdefault((a.GetSymbol(), r), a.GetIdx())

        # Charge just those representatives
        for atom_idx in reps.values():
            mol_copy = RWMol(mol)
            mol_copy.GetAtomWithIdx(atom_idx).SetFormalCharge(CHARGE)
            try:
                Chem.SanitizeMol(mol_copy)
                csmi = Chem.MolToSmiles(mol_copy, canonical=True)
            except Exception:
                continue

            if csmi not in seen:
                seen.add(csmi)
                charged.append(csmi)

    # Write unique, sorted
    path = os.path.join(OUTPUT_DIR, f"chargedSMILES_{FORMULA}.txt")
    with open(path, "w") as f:
        for c in sorted(charged):
            f.write(c + "\n")
    return charged


def filter_by_mass(smiles_list):
    filtered = []
    seen = set()
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if not mol:
            continue
        mass = Descriptors.ExactMolWt(mol)
        if abs(mass - TARGET_MASS) <= TOLERANCE:
            csmi = Chem.MolToSmiles(mol, canonical=True)
            if csmi not in seen:
                seen.add(csmi)
                filtered.append(csmi)

    path = os.path.join(OUTPUT_DIR, f"filteredchargedSMILES_{FORMULA}.txt")
    with open(path, "w") as f:
        for c in sorted(filtered):
            f.write(c + "\n")
    return filtered


######################### OUTPUT REPORTING ##########################################
def generate_xyz_from_smiles(smiles):
    """Generate optimized XYZ coordinates from SMILES string"""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"[WARNING] Could not parse SMILES: {smiles}")
            return None
            
        # Add hydrogens while preserving charges
        mol_with_h = Chem.AddHs(mol)
        
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

        mol_with_h = Chem.AddHs(mol)

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

        formula = CalcMolFormula(mol)
        mass = f"{Descriptors.ExactMolWt(mol):.4f}"
        
        # Always generate PNG with annotations
        try:
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

########################### SUMMARY FILE WRITER #################################
def _safe_count_lines(path):
    try:
        with open(path, "r") as f:
            return sum(1 for line in f if line.strip())
    except Exception:
        return 0

def _write_run_summary(gen_time_s, gen_ram_mb, charge_time_s, charge_ram_mb,
                       total_time_s, total_ram_mb):
    """
    Writes a human-readable summary file mirroring the example the user shared.
    Only uses information already written to disk by the normal pipeline.
    """
    nsmiles_file = f"nSMILES_{FORMULA}.txt"
    charged_file = f"chargedSMILES_{FORMULA}.txt"
    filtered_file = f"filteredchargedSMILES_{FORMULA}.txt"

    nsmiles_count = _safe_count_lines(os.path.join(OUTPUT_DIR, nsmiles_file))
    filtered_count = _safe_count_lines(os.path.join(OUTPUT_DIR, filtered_file))

    # Build summary text
    lines = []
    lines.append("ForMileS Run Summary")
    lines.append("====================\n")
    lines.append(f"Formula: {FORMULA}")
    lines.append(f"Charge: {CHARGE}")
    lines.append(f"Target mass: {TARGET_MASS} ± {TOLERANCE}\n")
    lines.append("Timings & Memory")
    lines.append("-----------------")
    lines.append(f"Molecular graph generation:   time = {gen_time_s:.3f} s | RAM ~ {gen_ram_mb:.1f} MB" if gen_ram_mb is not None else
                 f"Molecular graph generation:   time = {gen_time_s:.3f} s | RAM ~ unavailable")
    lines.append(f"Charge generation:            time = {charge_time_s:.3f} s | RAM ~ {charge_ram_mb:.1f} MB" if charge_ram_mb is not None else
                 f"Charge generation:            time = {charge_time_s:.3f} s | RAM ~ unavailable")
    lines.append(f"Total wall-time:              time = {total_time_s:.3f} s | RAM ~ {total_ram_mb:.1f} MB\n" if total_ram_mb is not None else
                 f"Total wall-time:              time = {total_time_s:.3f} s | RAM ~ unavailable\n")
    lines.append("Counts")
    lines.append("------")
    lines.append(f"nSMILES file count:                 {nsmiles_count}")
    lines.append(f"filteredchargedSMILES file count:   {filtered_count}\n")
    lines.append("Files")
    lines.append("-----")
    lines.append(f"nSMILES file:               {nsmiles_file}")
    lines.append(f"charged file:               {charged_file}")
    lines.append(f"filtered charged file:      {filtered_file}\n")
    lines.append("Notes")
    lines.append("-----")
    lines.append("RAM values are instantaneous RSS snapshots taken right after each step.\n")

    # File name mirrors example: run_summary_{FORMULA}.txt
    out_path = os.path.join(OUTPUT_DIR, f"run_summary_{FORMULA}.txt")
    with open(out_path, "w") as f:
        f.write("\n".join(lines))

    print(f"[OK] Wrote run summary: {out_path}")

########################### MAIN EXECUTION ######################################
if __name__ == "__main__":
    t0 = time.perf_counter()
    
    print("================== ForMileS v2.5 ==================")
    print(f"Structural Setup: Molecular Branching={ALLOW_BRANCHING}, Cyclic={ALLOW_CYCLES}")
    create_output_folder()

    # --- Stage 1: Molecular graph generation (timing + RAM snapshot)
    t_gen0 = time.perf_counter()
    base_smiles = run_generation()
    gen_time = time.perf_counter() - t_gen0
    gen_ram = _get_ram_mb()

    # --- Stage 2: Charge generation (timing + RAM snapshot)
    t_c0 = time.perf_counter()
    charged = generate_charged_smiles(base_smiles)
    charge_time = time.perf_counter() - t_c0
    charge_ram = _get_ram_mb()

    # --- Stage 3: Mass filter (no special timing requested by user)
    final = filter_by_mass(charged)

    # --- Images / coordinate outputs (unchanged behavior)
    smiles_to_images(final)

    print("=================== FINISHED ;) ===================")

    # NEW: end timer + RAM and write summary file
    elapsed = time.perf_counter() - t0
    ram_mb = _get_ram_mb()
    if ram_mb is not None:
        print(f"[SUMMARY] Elapsed={elapsed:.2f}s | RAM={ram_mb:.1f} MB (RSS)")
    else:
        print(f"[SUMMARY] Elapsed={elapsed:.2f}s | RAM=unavailable (install 'psutil' to enable)")

    # Persist summary mirroring user's example
    _write_run_summary(
        gen_time_s=gen_time,
        gen_ram_mb=(gen_ram if gen_ram is not None else float("nan")),
        charge_time_s=charge_time,
        charge_ram_mb=(charge_ram if charge_ram is not None else float("nan")),
        total_time_s=elapsed,
        total_ram_mb=(ram_mb if ram_mb is not None else float("nan")),
    )

