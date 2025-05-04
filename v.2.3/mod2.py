#
# This module is designed to filter the combinatory space of SMILES
# generated in mod1.py. 
# The filter must to be prepared as such as to incorporte the most
# important features of the precursor molecule.
#
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdqueries
from rdkit.Chem.rdMolDescriptors import CalcNumRings
from tqdm import tqdm
import os


def is_linear_chain(mol):
    """Check if molecule is a linear (non-cyclic, non-branched) chain."""
    # Must be a tree (no cycles)
    if mol.GetRingInfo().NumRings() > 0:
        return False

    # Linear means: only two atoms with degree 1 (ends), all others degree 2
    degrees = [atom.GetDegree() for atom in mol.GetAtoms()]
    num_ends = degrees.count(1)
    num_middle = degrees.count(2)

    return num_ends == 2 and (num_ends + num_middle == len(degrees))


def contains_COC_bond(mol):
    """Check for C-O-C using a SMARTS pattern."""
    coc_pattern = Chem.MolFromSmarts("C-C-O-C-C")
    return mol.HasSubstructMatch(coc_pattern)


def filter_smiles(input_file, FORMULA, output_file=None):

    output_dir = f"OutputFiles_{FORMULA}"
    output_path = os.path.join(output_dir, output_file)
    input_path = os.path.join(output_dir, input_file)

    if output_file is None:
        output_file = f"ParentRelatedSMILES_{FORMULA}.txt"
    with open(input_path, "r") as f:
        smiles_list = [line.strip() for line in f if line.strip()]

    filtered = []

    for smi in tqdm(smiles_list, desc="Filtering SMILES"):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        if contains_COC_bond(mol) and is_linear_chain(mol):
            filtered.append(smi)

    # Save results
    with open(output_path, "w") as f:
        for smi in filtered:
            f.write(smi + "\n")

    print(f"Saved {len(filtered)} filtered SMILES to '{output_path}'")