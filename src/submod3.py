#########################################################################
# This module is designed to generate all possible charged SMILES for a #
# set of neutral SMILES that the user already filtered according to the #
# features of importance in the precursor. Also, a last important       #
# filter of mass garantee the exact mass of the precursor ion           #
#########################################################################


import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import math
import os

from rdkit import Chem
from rdkit.Chem import Descriptors
from itertools import combinations, product
from rdkit.Chem import Draw
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from tqdm import tqdm 

########## generate charged SMILES from ParentRelatedSMILES ##########

def generate_charged_smiles(FORMULA, CHARGE, input_file, output_file=None):

    output_dir = f"OutputFiles_{FORMULA}_Charge_{CHARGE}" 
    output_path = os.path.join(output_dir, output_file)
    input_path = os.path.join(output_dir, input_file)

    if output_file is None:
        output_file = f"chargedSMILES_{FORMULA}.txt"
    with open(input_path, "r") as f:
        smiles_list = [line.strip() for line in f if line.strip()]

    charged_smiles = []

    for smi in tqdm(smiles_list, desc="Building charged SMILES"):
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

    with open(output_path, "w") as f:
        for smi in charged_smiles:
           f.write(smi + "\n")

    print(f"Saved {len(charged_smiles)} charged SMILES to '{output_path}'")

################ Filter the charged SMILES by mass ################

def filter_charged_smiles_by_mass(FORMULA, CHARGE, input_file, TARGET_MASS, tolerance, output_file=None):
    
    output_dir = f"OutputFiles_{FORMULA}_Charge_{CHARGE}" 
    output_path = os.path.join(output_dir, output_file)
    input_path = os.path.join(output_dir, input_file)
    
    if output_file is None:
        output_file = f"filteredchargedSMILES_{FORMULA}.txt"
    with open(input_path, "r") as f:
        smiles_list = [line.strip() for line in f if line.strip()]

    filtered_charged = []
    
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            mass = Descriptors.ExactMolWt(mol)
            if abs(mass - TARGET_MASS) <= tolerance:
                filtered_charged.append(smi)

    with open(output_path, "w") as f:
        for smi in filtered_charged:
            f.write(smi + "\n")
    print(f"Saved {len(filtered_charged)} filtered charged SMILES to '{output_path}'")