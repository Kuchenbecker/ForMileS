# ForMileS v2.5. 
 -- acronymn for Formation of Mass SMILES --

Created by Vinicius Kuchenbecker (1) in association with UNICAMP.
(1) vincius.kuchenbecker@gmail.com

In ESI-Tandem MS/MS experiments, one may know the precursor molecular formula, charge state (generally positive) and structure.
After collision, MS spectra will give two main informations: Exact mass of fragment ion and charge state (positive or negative).

This program is designed to survey the combinatory space of isomers generated from this two informations:
- **1.Proportion of heavy atoms in the fragment molecule. That is, the molecular formula withou explicit hydrogens**
- **2.A estructural chemical feature of the precurson ion that must be in the fragment ion as well**
- **3.The charge state of fragment ion**
- **4.The Exact mass of fragment ion**    
Hence, it is suited for cases in which one is in position of parsing this informations.

#############################################################################

Inputs of the program:
1. FORMULA
2. PRECURSOR_FEATURE
3. CHARGE
4. EXACT_MASS

User must to give this inside **main.py** script.

**main.py** will work calling four separated modules that will generate the final wanted result.

For futher detailing how each submodule works and what should or not be altered, refer to documentation.
Some comments about usage are also within each module or submodule.

This version is configured to be used with C and O positive charged molecules. Also is configured in Module 1 and 2 to avoid cyclic molecules. Don't forget to check it.

Any comment, help or suggestion, adress the e-mail.

