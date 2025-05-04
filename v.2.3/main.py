# ----- LIST OF MODULES AND FUNCTIONS IMPORTED
import os
from mod1 import generate_smiles
from mod2 import filter_smiles
from mod3 import generate_charged_smiles
from mod3 import filter_charged_smiles_by_mass
from mod4 import smiles_to_molecules

# ----- Inputs

FORMULA = "C6O2"     #molecular formula without explicit hydrogens
TARGET_MASS = 117.092     #Charged fragment target mass

# ----- output folder creation
def create_output_folder(FORMULA):  
    script_dir = os.path.dirname(os.path.abspath(__file__))
    folder_name = f"OutputFiles_{FORMULA}"
    folder_path = os.path.join(script_dir, folder_name)

    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
        print(f"Folder created: {folder_path}")
    else:
        print(f"Folder already exists: {folder_path}")
    
    return folder_path

create_output_folder(FORMULA)

# ----- Main functions and ouputs/prints
print("###############################################################")
print("#################   ForMileS v.2.3 Starting   #################")
print("###############################################################")
print("                                                               ")
print("                                                               ")

print(f"Generating neutral SMILES for: {FORMULA}")
generate_smiles(FORMULA, f"nSMILES_{FORMULA}.txt")

print("Checking heritage with precursor molecule")
print("Generating filtered by heritage SMILES")
filter_smiles(f"nSMILES_{FORMULA}.txt", FORMULA, f"ParentRelatedSMILES_{FORMULA}.txt")
 
print("Generating charged SMILES from ParentRelatedSMILES")
print("Filtering charged SMILES by mass")
generate_charged_smiles(FORMULA, f"ParentRelatedSMILES_{FORMULA}.txt", f"chargedSMILES_{FORMULA}.txt")
filter_charged_smiles_by_mass(FORMULA, f"chargedSMILES_{FORMULA}.txt", TARGET_MASS, 0.05, f"filteredchargedSMILES_{FORMULA}.txt")

print("Generating charged molecules properties calculations and images")
smiles_to_molecules(FORMULA, f"filteredchargedSMILES_{FORMULA}.txt")
print("All images saved to OutputFiles folder")

print("                                                               ")
print("                                                               ")
print("###############################################################")
print("#############   ForMileS v.2.3 Ended Succesfully   ############")
print("###############################################################")
