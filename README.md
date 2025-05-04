# ForMileS v.2.3 
 -- Formation of Mass Smiles --

Created by Vinicius Kuchenbecker (1) in association with UNICAMP.
(1) vincius.kuchenbecker@gmail.com

In ESI-Tandem MS/MS experiments, one may know the precursor molecular formula, charge state (generally positive) and structure.
After collision, MS spectra will give two main informations: Exact mass of fragment ion and charge state (positive or negative).

This program is designed to survey the combinatory space of isomers generated from this two informations:
- **1.Exact mass of charged fragment**
- **2.(proposed) molecular formula of fragment ion**
    
Hence, it is suited for cases in which one is in position of parsing this two input informations.

#################################################################################

Inputs of the program:
1. FORMULA
2. EXACT_MASS

User must to give this inside **main.py** script.

**main.py** will work calling four separated modules that will generate the final wanted result.

# Module 1 
Module 1 generates a extense list of possible structures as SMILES based on the molecular formula (without hydrogens for simplicity of script).
It will do so based on simple valence rules that can be configured in Module 1. Graph theory is the main core of module 1, in which valence is tested assuming atoms are nodes and bonds are arest.
One may be prone to configure this module to consider every possible cyclic molecule, or not. All possible bond orders, or not... and so on.

In this version, the focus is C and O molecules only, without cyclic generation. Bear in mind that, the bigger the molecule and the less restrictions, the longer the calculation of all possibilities.

# Module 2
Module 2 filters the Module 1. Specifically, this Module will look up to features to stablish a heritage relationship between precursor molecule knowlodge and fragment.
The filter is configured the way the user wants. It is up to the user to say what is considered an essential and non-essential chemical feature of precursor ion that must be in fragment ion.

The way the Module 2 understand this is using SMARTS rules.

In simpler words: In this step, a large set of combinatory structures become a smaller set of molecules that are possible "children" of the precursor.

# Module 3
As stated in the introduction, ESI-MS/MS will always involve both precursor and fragment charged.

The Module 3 will:
1º: Generate all possible charged isomers in every acceptle atomic position (respecting some basic valence rule) from the list generated in Module 2.
2º: Filter only the charged isomers that have the exact mass parsed by the user as input.

Here is where the second experimental MS spectra information acts as filter too. 

# Module 4

Just for printing the final structures in 2D with exact mass. 

## Some other comments about the execution

1. **main.py** call all the modules to work in the right sequence.
2. ouput files stay always in a folder generated when **main.py** is called.
3. One can alter **main.py** to, for exemple, skip execution of one module. But be aware of the interdepencies.
4. In the Output Folder, SMILES from each module execution are saved as .txt files. Images as .jpeg are generated only for the final list of charged SMILES.

This version is configured to be used with C and O positive charged molecules. Also is configured in Module 1 and 2 to avoid cyclic molecules. Don't forget to check it.

Any comment, help or suggestion, adress the e-mail.

