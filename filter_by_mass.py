from rdkit import Chem
from rdkit.Chem import Descriptors
from smiles_generator import generate_smiles  # assumes the first script is named smiles_generator.py

TARGET_MASS = 60.058
TOLERANCE = 0.0005  # Allows for floating-point rounding tolerance

def filter_smiles_by_mass(smiles_list, target_mass):
    filtered = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            mass = Descriptors.ExactMolWt(mol)
            if abs(mass - target_mass) <= TOLERANCE:
                filtered.append(smi)
    return filtered

if __name__ == "__main__":
    formula = "C6O3"
    smiles_list = generate_smiles(formula)
    filtered_smiles = filter_smiles_by_mass(smiles_list, TARGET_MASS)
    
    print("Filtered SMILES (Exact Mass ~60.058):")
    for smi in filtered_smiles:
        print(smi)

