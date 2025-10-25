# 🧪 ForMileS — Formation of Mass SMILES

**ForMileS (Formation of Mass SMILES)** is a molecular graph generation engine built for **systematic SMILES enumeration**, **mass-based filtering**, and **charge modeling**.  
It allows researchers to explore chemical space consistent with a given **molecular formula**, **exact mass**, and **bonding rules** using RDKit, while respecting valence and structural constraints.

<p align="center">
  <img src="https://github.com/Kuchenbecker/ForMileS/assets/banner_example.png" width="700"/>
</p>

---

## 🔬 Key Features

- **Automated Molecular Graph Expansion**  
  Expands molecular scaffolds recursively following valence rules, bond orders, and formula constraints.

- **Charge Enumeration**  
  Assigns formal charges to symmetry-unique atoms according to allowed elements.

- **Mass Filtering**  
  Retains only molecules within a specified tolerance of the target mass.

- **Structure Visualization and Export**
  - PNG and optional SVG 2D structure images with formula/mass annotation.
  - 3D coordinate generation (`.xyz` and `.mol`) using RDKit + UFF optimization.

- **Parameter File (`parameters.json`)**
  - Defines bond orders, max valence per element, and allowed chargeable atoms.
  - Easy to customize for different chemical systems.

- **Optional GUI (`gui.py`)**
  - Simple cross-platform Tkinter interface for defining formula, scaffold, charge, and export options.
  - Allows structure generation without command-line usage.

---

## ⚙️ Installation

### Requirements

| Library | Purpose |
|----------|----------|
| `rdkit` | Molecule construction, SMILES, and 3D optimization |
| `pillow` | Image annotation and export |
| `tqdm` | Progress visualization |
| `psutil` *(optional)* | Memory usage monitoring |
| `tkinter` | GUI interface (usually included with Python) |

Install dependencies:

```bash
pip install rdkit-pypi pillow tqdm psutil
