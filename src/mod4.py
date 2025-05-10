from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors
import matplotlib.pyplot as plt
import os
import math

####### Image Generation Parameters #######
MOLS_PER_IMAGE = 10  # Max molecules per image file
MOLS_PER_ROW = 5     # Molecules per row in each image
IMAGE_SIZE = (40, 50)  # Base figure size (scaled per image)

def smiles_to_molecules(FORMULA, input_file):
    output_dir = f"OutputFiles_{FORMULA}"
    input_path = os.path.join(output_dir, input_file)
    
    with open(input_path, "r") as f:
        smiles_list = [line.strip() for line in f if line.strip()]

    molecules = []
    data = []

    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            mass = f"{Descriptors.ExactMolWt(mol):.4f}"
            formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
            data.append((smiles, mass, formula))
            molecules.append(mol)
        else:
            print(f"Warning: Could not parse SMILES '{smiles}'")

    if not molecules:
        print("No valid molecules to display.")
        return molecules, data

    # Split molecules into chunks of 10
    total_molecules = len(molecules)
    num_images = math.ceil(total_molecules / MOLS_PER_IMAGE)

    for img_idx in range(num_images):
        start_idx = img_idx * MOLS_PER_IMAGE
        end_idx = min((img_idx + 1) * MOLS_PER_IMAGE, total_molecules)
        chunk_mols = molecules[start_idx:end_idx]
        chunk_data = data[start_idx:end_idx]

        n_rows = math.ceil(len(chunk_mols) / MOLS_PER_ROW)
        fig, axes = plt.subplots(n_rows, MOLS_PER_ROW, 
                                figsize=(IMAGE_SIZE[0], IMAGE_SIZE[1] * (n_rows / 5)))
        
        if n_rows == 1:
            axes = [axes] if isinstance(axes, plt.Axes) else axes.reshape(1, -1)
        else:
            axes = axes.ravel()

        for ax in axes:
            ax.axis('off')

        for i, (mol, (smiles, mass, formula)) in enumerate(zip(chunk_mols, chunk_data)):
            if i >= len(axes):
                break
                
            img = Draw.MolToImage(mol, size=(250, 100))
            axes[i].imshow(img)
            axes[i].axis('off')
            
            annotation = f"{smiles}\nMass: {mass}\n{formula}"
            axes[i].text(0.5, 0.01, annotation,
                        transform=axes[i].transAxes,
                        ha='center', va='bottom',
                        fontsize=6,
                        bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))

        plt.tight_layout()
        output_path = os.path.join(output_dir, f"MolFormile_Plot_{FORMULA}_{img_idx + 1}.jpeg")
        plt.savefig(output_path, format='jpeg', dpi=300)
        plt.close()
        print(f"Saved image {img_idx + 1}/{num_images} to '{output_path}'")

    return molecules, data