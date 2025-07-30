#########################################################################
# This module is designed to filter the combinatory space of SMILES     #
# generated in mod1.py.                                                 #
# The filter must to be prepared as such as to incorporte the most      #
# important features of the precursor molecule.                         #
#########################################################################

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from tqdm import tqdm
import os


def is_linear(mol):
    """
    Check if molecule is a linear (non-cyclic, non-branched) chain.
    Conditions:
    - Acyclic (no rings)
    - Exactly two terminal atoms (degree 1), all others degree 2
    """
    degrees = [atom.GetDegree() for atom in mol.GetAtoms()]
    num_ends = degrees.count(1)
    num_middle = degrees.count(2)

    return num_ends == 2 and (num_ends + num_middle == len(degrees))


def contains_feature(mol, PRECURSOR_FEATURES):
    """Check for precursor chemical features using SMARTS rules."""
    for feature in PRECURSOR_FEATURES:
        pattern = Chem.MolFromSmarts(feature)
        if pattern and mol.HasSubstructMatch(pattern):
            return True
    return False


def filter_smiles(input_file, FORMULA, CHARGE, PRECURSOR_FEATURES, branched, ring, output_file=None):
    """
    Filters SMILES strings based on precursor features and structure rules.

    Parameters:
    - branched (bool): If False, only linear (non-branched) chains are accepted.
    - ring (bool): If False, cyclic molecules are rejected.
    """
    output_dir = f"OutputFiles_{FORMULA}_Charge_{CHARGE}"
    output_path = os.path.join(output_dir, output_file if output_file else f"ParentRelatedSMILES_{FORMULA}.txt")
    input_path = os.path.join(output_dir, input_file)

    with open(input_path, "r") as f:
        smiles_list = [line.strip() for line in f if line.strip()]

    filtered = []

    for smi in tqdm(smiles_list, desc="Filtering SMILES"):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue

        if not ring and mol.GetRingInfo().NumRings() > 0:
            continue

        if not branched and not is_linear(mol):
            continue

        if contains_feature(mol, PRECURSOR_FEATURES):
            filtered.append(smi)

    with open(output_path, "w") as f:
        for smi in filtered:
            f.write(smi + "\n")

    print(f"Saved {len(filtered)} filtered SMILES to '{output_path}'")