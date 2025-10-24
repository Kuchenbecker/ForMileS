
#################################################################################
###                   ForMileS v3.1 (Explicit-H Heavy-First, fixes)          ###
#################################################################################

import os, sys, json, time
from collections import Counter
from rdkit import Chem
from rdkit.Chem.rdchem import BondType

CONFIG_FILE = "config.json"
PARAMS_FILE = "parameters.json"

def load_json(path, default=None):
    try:
        with open(path, "r") as f:
            return json.load(f)
    except Exception:
        return default

def resource_path(relative_path):
    try:
        base_path = sys._MEIPASS  # type: ignore
    except Exception:
        base_path = os.path.abspath(".")
    return os.path.join(base_path, relative_path)

def parse_formula(formula_str: str) -> Counter:
    import re
    pat = re.compile(r"([A-Z][a-z]?)(\d*)")
    c = Counter()
    for elem, num in pat.findall(formula_str):
        c[elem] += int(num) if num else 1
    return c

def count_atoms(mol: Chem.Mol) -> Counter:
    c = Counter()
    for a in mol.GetAtoms():
        c[a.GetSymbol()] += 1
    return c

def norm_bokey(a: str, b: str) -> str:
    """Make key like the JSON uses: ["A","B"] with double quotes and comma, sorted."""
    a, b = sorted([a, b])
    return f'["{a}","{b}"]'

def get_allowed_bond_orders(params, sym1, sym2):
    key = norm_bokey(sym1, sym2)
    return params.get("bond_orders", {}).get(key, [])

def open_sites(mol: Chem.Mol, max_valence: dict):
    sites = []
    for a in mol.GetAtoms():
        sym = a.GetSymbol()
        vmax = max_valence.get(sym, 4)
        # current explicit valence (sum bond orders)
        val = 0
        for b in a.GetBonds():
            val += int(b.GetBondTypeAsDouble())
        # Consider that implicit H may still fill up to vmax
        if val < vmax:
            sites.append(a.GetIdx())
    return sites

def attach_atom(mol: Chem.Mol, sym: str, anchor_idx: int, order: int):
    rw = Chem.RWMol(mol)
    new_idx = rw.AddAtom(Chem.Atom(sym))
    rw.AddBond(anchor_idx, new_idx, [BondType.SINGLE, BondType.DOUBLE, BondType.TRIPLE][order-1])
    nm = rw.GetMol()
    # sanitize properties only (avoid over-restrictive kekulization mid-growth)
    Chem.SanitizeMol(nm, sanitizeOps=Chem.SanitizeFlags.SANITIZE_PROPERTIES)
    return nm

def count_implicit_H_needed(mol: Chem.Mol) -> int:
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        pass
    total = 0
    for a in mol.GetAtoms():
        h = a.GetNumImplicitHs()
        if h < 0: h = 0
        total += h
    return total

def decorate_with_exact_H(mol: Chem.Mol, target_H: int):
    m = Chem.Mol(mol)
    try:
        mH = Chem.AddHs(m, addCoords=False)
        h_count = sum(1 for a in mH.GetAtoms() if a.GetSymbol() == "H")
        if h_count != target_H:
            return None
        return mH
    except Exception:
        return None

def smiles_with_explicit_H(m: Chem.Mol) -> str:
    return Chem.MolToSmiles(m, canonical=True, allHsExplicit=True, allBondsExplicit=False)

DEFAULTS = {
    "ALLOW_CYCLES": True,
    "ALLOW_BRANCHING": True,
    "MAX_DOUBLE_BONDS": 12,
    "MAX_TRIPLE_BONDS": 6,
    "OUTPUT_DIR": "output",
    # keys expected from config: FORMULA, MOL_SCAFFOLD (list of SMILES), OUTPUT_DIR, etc.
}

class ForMileSGenerator:
    def __init__(self, config, params):
        self.cfg = {**DEFAULTS, **(config or {})}
        formula = self.cfg.get("FORMULA") or self.cfg.get("TARGET_FORMULA") or "C3H8O"
        self.target = parse_formula(formula)
        self.target_H = self.target.get("H", 0)
        self.target_nonH = Counter({k:v for k,v in self.target.items() if k != "H"})
        self.params = params or {}
        self.max_valence = self.params.get("max_valence", {})
        self.collected, self.seen = [], set()
        outdir = self.cfg.get("OUTPUT_DIR", "output")
        os.makedirs(outdir, exist_ok=True)
        self.outdir = outdir

    def deficit(self, mol: Chem.Mol) -> Counter:
        current = count_atoms(mol)
        d = Counter()
        for k,v in self.target_nonH.items():
            d[k] = max(0, v - current.get(k,0))
        return d

    def heavy_complete(self, mol: Chem.Mol) -> bool:
        c = count_atoms(mol)
        for k,v in self.target_nonH.items():
            if c.get(k,0) != v: return False
        return True

    def grow(self, mol: Chem.Mol):
        # Bound by minimal implicit H
        minH = count_implicit_H_needed(mol)
        if minH > self.target_H:
            return

        if self.heavy_complete(mol):
            self.try_accept(mol)
            return

        d = self.deficit(mol)
        if not d: return

        # try rarer elements first
        add_list = sorted([k for k in d.keys() if d[k] > 0], key=lambda x: d[x])
        sites = open_sites(mol, self.max_valence)
        if not sites: return

        progressed = False
        for sym in add_list:
            for idx in sites:
                anchor_sym = mol.GetAtomWithIdx(idx).GetSymbol()
                allowed = get_allowed_bond_orders(self.params, anchor_sym, sym)
                if not allowed:
                    continue
                for order in allowed:
                    try:
                        new_m = attach_atom(mol, sym, idx, order)
                        self.grow(new_m)
                        progressed = True
                        if len(self.collected) >= self.cfg.get("MAX_STRUCTURES", 10000):
                            return
                    except Exception:
                        continue
        if not progressed:
            return

    def try_accept(self, mol: Chem.Mol):
        mH = decorate_with_exact_H(mol, self.target_H)
        if mH is None: return
        smi = smiles_with_explicit_H(mH)
        if not smi or smi in self.seen: return
        self.seen.add(smi)
        self.collected.append(smi)

    def seed_molecule(self) -> Chem.Mol:
        scaffolds = self.cfg.get("MOL_SCAFFOLD") or []
        if isinstance(scaffolds, list) and scaffolds:
            # take first scaffold SMILES and keep only heavy atoms (H stripped)
            m = Chem.MolFromSmiles(scaffolds[0])
            if m is not None:
                m = Chem.RemoveHs(m)
                # prune atoms not in target_nonH
                rw = Chem.RWMol()
                map_old_to_new = {}
                for a in m.GetAtoms():
                    if a.GetSymbol() in self.target_nonH:
                        idx = rw.AddAtom(Chem.Atom(a.GetSymbol()))
                        map_old_to_new[a.GetIdx()] = idx
                # add bonds when both ends kept
                for b in m.GetBonds():
                    i = b.GetBeginAtomIdx(); j = b.GetEndAtomIdx()
                    if i in map_old_to_new and j in map_old_to_new:
                        rw.AddBond(map_old_to_new[i], map_old_to_new[j], b.GetBondType())
                try:
                    nm = rw.GetMol()
                    Chem.SanitizeMol(nm, sanitizeOps=Chem.SanitizeFlags.SANITIZE_PROPERTIES)
                    return nm
                except Exception:
                    pass
        # fallback: start with any heavy element from target
        if self.target_nonH:
            sym = next(iter(self.target_nonH.keys()))
            rw = Chem.RWMol(); rw.AddAtom(Chem.Atom(sym))
            return rw.GetMol()
        # hydrogen-only (rare)
        return Chem.MolFromSmiles("[H]")

    def run(self):
        seed = self.seed_molecule()
        self.grow(seed)
        out_smiles = os.path.join(self.outdir, "nSMILES_explicitH.txt")
        with open(out_smiles, "w") as f:
            for s in self.collected:
                f.write(s + "\n")
        print(f"[ForMileS Explicit-H] Generated {len(self.collected)} structures -> {out_smiles}")
        return out_smiles, len(self.collected)

def main():
    if getattr(sys, "frozen", False):
        os.chdir(os.path.dirname(sys.executable))
    else:
        os.chdir(os.path.dirname(os.path.abspath(__file__)))

    cfg = load_json(resource_path(CONFIG_FILE), {})
    params = load_json(resource_path(PARAMS_FILE), {})
    gen = ForMileSGenerator(cfg, params)
    gen.run()

if __name__ == "__main__":
    main()
