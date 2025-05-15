#################################################################################
###                                                                           ###
###                   ForMileS - FORMATION OF MASS SMILES                     ###
###                                                                           ###
#################################################################################
# If the user wants to disable one of the submodules, they can comment out      #
# the import statement and the corresponding function call in the main module.  #
# This will prevent the submodule from being executed, but the rest of the code #
# will still run.                                                               #
#################################################################################

########################## IMPORT MODULES OF FORMILES ###########################
import os
from submod1 import generate_smiles
from submod2 import filter_smiles
from submod3 import generate_charged_smiles
from submod3 import filter_charged_smiles_by_mass
from submod4 import smiles_to_molecules

############################### INPUT PARAMETERS ################################
#
#
#
FORMULA = "C4O2"
PRECURSOR_FEATURES = ["C-O-C-C-C", "C-O-C-C-C-O"]
CHARGE = +1
TARGET_MASS = 89.060
#
#
#
############################# DIRECTORY CREATION ###############################
def create_output_folder(FORMULA, CHARGE):  
    script_dir = os.path.dirname(os.path.abspath(__file__))
    folder_name = f"OutputFiles_{FORMULA}_Charge_{CHARGE}"
    folder_path = os.path.join(script_dir, folder_name)

    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
        print(f"Folder created: {folder_path}")
    else:
        print(f"Folder already exists: {folder_path}")
    
    return folder_path

create_output_folder(FORMULA, CHARGE)

###################### EXECUTION OF FUNCTIONS PER MODULE #######################
print("###############################################################")
print("#################   ForMileS v.2.3 Starting   #################")
print("###############################################################")
print("                                                               ")
print("                                                               ")

print(f"Generating neutral SMILES for: {FORMULA}")
generate_smiles(FORMULA, CHARGE, f"nSMILES_{FORMULA}.txt")

print("Checking heritage with precursor molecule")
print("Generating filtered by heritage SMILES")
filter_smiles(f"nSMILES_{FORMULA}.txt", FORMULA, CHARGE, PRECURSOR_FEATURES, f"ParentRelatedSMILES_{FORMULA}.txt")
 
print("Generating charged SMILES from ParentRelatedSMILES")
print("Filtering charged SMILES by mass")
generate_charged_smiles(FORMULA, CHARGE, f"ParentRelatedSMILES_{FORMULA}.txt", f"chargedSMILES_{FORMULA}.txt")
filter_charged_smiles_by_mass(FORMULA, CHARGE, f"chargedSMILES_{FORMULA}.txt", TARGET_MASS, 0.05, f"filteredchargedSMILES_{FORMULA}.txt")

print("Generating charged molecules properties calculations and images")
smiles_to_molecules(FORMULA, CHARGE, f"filteredchargedSMILES_{FORMULA}.txt")
print("All images saved to OutputFiles folder")

print("                                                               ")
print("                                                               ")
print("###############################################################")
print("#############   ForMileS v.2.3 Ended Succesfully   ############")
print("###############################################################")
