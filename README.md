# ForMileS
-- acronymn for Formation of Mass SMILES --

Author: Vinicius Kuchenbecker

In ESI-Tandem MS/MS experiments, one may know the precursor molecular formula, charge state (generally positive) and structure.
After collision, MS spectra will give two main informations: Exact mass of fragment ion and charge state (positive or negative).

This program is designed to survey the combinatory space of isomers generated from this two informations:
- **1.Proportion of heavy atoms in the fragment molecule. That is, the molecular formula withou explicit hydrogens**
- **2.A estructural chemical feature of the precurson ion that must be in the fragment ion as well**
- **3.The charge state of fragment ion**
- **4.The Exact mass of fragment ion**    
Hence, it is suited for cases in which one is in position of parsing this informations.

#############################################################################

Version 1: Uses extensive numeration method - No GUI

Version 2: Branch-and-Bound (B&B) with SMARTS - No GUI (CURRENT ONE)

Version 3: B&B with GUI (IN DEV)
