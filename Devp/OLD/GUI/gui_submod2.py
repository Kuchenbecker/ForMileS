from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import os

def is_linear(mol):
    degrees = [atom.GetDegree() for atom in mol.GetAtoms()]
    num_ends = degrees.count(1)
    num_middle = degrees.count(2)
    return num_ends == 2 and (num_ends + num_middle == len(degrees))

def contains_feature(mol, PRECURSOR_FEATURES):
    for feature in PRECURSOR_FEATURES:
        pattern = Chem.MolFromSmarts(feature)
        if pattern and mol.HasSubstructMatch(pattern):
            return True
    return False

def filter_smiles(input_file, FORMULA, CHARGE, PRECURSOR_FEATURES, branched, ring, output_file=None,
                 should_stop=None, update_callback=None):
    output_dir = f"OutputFiles_{FORMULA}_Charge_{CHARGE}"
    output_path = os.path.join(output_dir, output_file or f"ParentRelatedSMILES_{FORMULA}.txt")
    input_path = os.path.join(output_dir, input_file)

    with open(input_path, "r") as f:
        smiles_list = [line.strip() for line in f if line.strip()]

    filtered = []
    total = len(smiles_list)
    BATCH_SIZE = 1000
    UPDATE_FREQ = max(100, total // 100)

    for i in range(0, total, BATCH_SIZE):
        batch = smiles_list[i:i+BATCH_SIZE]
        for smi in batch:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue

            if not ring and mol.GetRingInfo().NumRings() > 0:
                continue
            if not branched and not is_linear(mol):
                continue
            if contains_feature(mol, PRECURSOR_FEATURES):
                filtered.append(smi)

        if should_stop and should_stop():
            break
            
        if update_callback:
            update_callback(min(i + BATCH_SIZE, total), total)

    with open(output_path, "w") as f:
        f.write("\n".join(filtered))
    
    print(f"Saved {len(filtered)} filtered SMILES to '{output_path}'")