from rdkit import Chem
from rdkit.Chem import Descriptors
import os

def generate_charged_smiles(FORMULA, CHARGE, input_file, output_file=None,
                            should_stop=lambda: False, update_callback=None):
    output_dir = f"OutputFiles_{FORMULA}_Charge_{CHARGE}"
    output_path = os.path.join(output_dir, output_file or f"chargedSMILES_{FORMULA}.txt")
    input_path = os.path.join(output_dir, input_file)

    with open(input_path, "r") as f:
        smiles_list = [line.strip() for line in f if line.strip()]

    charged_smiles = []
    total = len(smiles_list)

    for i, smi in enumerate(smiles_list):
        if should_stop():
            print("Stopped during charged SMILES.")
            break

        mol = Chem.MolFromSmiles(smi, sanitize=False)
        if mol is None:
            continue
        for atom in mol.GetAtoms():
            if atom.GetSymbol() in ["C", "O"]:
                mol_copy = Chem.RWMol(mol)
                atom_copy = mol_copy.GetAtomWithIdx(atom.GetIdx())
                atom_copy.SetFormalCharge(CHARGE)
                try:
                    charged_smi = Chem.MolToSmiles(mol_copy, canonical=True)
                    charged_smiles.append(charged_smi)
                except:
                    continue

        if update_callback and (i % 10 == 0 or i == total - 1):
            update_callback(i + 1, total)

    with open(output_path, "w") as f:
        for smi in charged_smiles:
            f.write(smi + "\n")
    print(f"Saved {len(charged_smiles)} charged SMILES to '{output_path}'")

def filter_charged_smiles_by_mass(FORMULA, CHARGE, input_file, TARGET_MASS, tolerance, output_file=None,
                                   should_stop=lambda: False, update_callback=None):
    output_dir = f"OutputFiles_{FORMULA}_Charge_{CHARGE}"
    output_path = os.path.join(output_dir, output_file or f"filteredchargedSMILES_{FORMULA}.txt")
    input_path = os.path.join(output_dir, input_file)

    with open(input_path, "r") as f:
        smiles_list = [line.strip() for line in f if line.strip()]

    filtered = []
    total = len(smiles_list)

    for i, smi in enumerate(smiles_list):
        if should_stop():
            print("Stopped during mass filtering.")
            break

        mol = Chem.MolFromSmiles(smi)
        if mol:
            mass = Descriptors.ExactMolWt(mol)
            if abs(mass - TARGET_MASS) <= tolerance:
                filtered.append(smi)

        if update_callback and (i % 10 == 0 or i == total - 1):
            update_callback(i + 1, total)

    with open(output_path, "w") as f:
        for smi in filtered:
            f.write(smi + "\n")
    print(f"Saved {len(filtered)} filtered charged SMILES to '{output_path}'")