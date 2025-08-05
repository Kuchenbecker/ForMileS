
import tkinter as tk
from tkinter import ttk, messagebox
import json
import subprocess
import os
import sys

CONFIG_FILE = "config.json"

def get_output_folder():
    try:
        with open(CONFIG_FILE) as f:
            config = json.load(f)
        formula = config.get("FORMULA", "Output")
        charge = config.get("CHARGE", 0)
        return f"OutputFiles_{formula}_Charge_{charge}"
    except:
        return "."

def run_program():
    try:
        config = {
            "FORMULA": formula_var.get(),
            "MOL_SCAFFOLD": [scaffold_var.get()],
            "CHARGE": int(charge_var.get()),
            "TARGET_MASS": float(mass_var.get()),
            "TOLERANCE": float(tolerance_var.get()),
            "PARAM_FILE": "parameters.json",
            "SAVE_XYZ": save_xyz_var.get(),
            "SAVE_MOL": save_mol_var.get(),
            "SAVE_SVG": save_svg_var.get(),
            "ALLOW_BRANCHING": branching_var.get(),
            "ALLOW_CYCLES": cycles_var.get(),
            "ALLOW_DOUBLE_BONDS": double_var.get(),
            "ALLOW_TRIPLE_BONDS": triple_var.get(),
            "MAX_DOUBLE_BONDS": int(max_double_var.get()),
            "MAX_TRIPLE_BONDS": int(max_triple_var.get())
        }

        with open(CONFIG_FILE, "w") as f:
            json.dump(config, f, indent=2)

        subprocess.run(["python", "main.py"])
        messagebox.showinfo("Finished", "ForMileS Ended Succefully. Check output directory.")
    except Exception as e:
        messagebox.showerror("Error", f"An error occurred trying to execute: {e}")

def open_output_folder():
    folder = get_output_folder()
    if not os.path.exists(folder):
        messagebox.showwarning("Directory not found", f"The directory folder {folder} was not created yet.")
        return

    if sys.platform == "win32":
        os.startfile(folder)
    elif sys.platform == "darwin":
        subprocess.Popen(["open", folder])
    else:
        subprocess.Popen(["xdg-open", folder])

root = tk.Tk()
root.title("ForMileS GUI")
root.geometry("400x630")

tk.Label(root, text="Molecular Formula:").pack()
formula_var = tk.StringVar(value="C7Cl")
tk.Entry(root, textvariable=formula_var).pack()

tk.Label(root, text="Basic Molecular Scaffold (BMS):").pack()
scaffold_var = tk.StringVar(value="C1C=CC=CC=1")
tk.Entry(root, textvariable=scaffold_var).pack()

tk.Label(root, text="Charge:").pack()
charge_var = tk.StringVar(value="-1")
tk.Entry(root, textvariable=charge_var).pack()

tk.Label(root, text="Exact Mass:").pack()
mass_var = tk.StringVar(value="125.016")
tk.Entry(root, textvariable=mass_var).pack()

tk.Label(root, text="Mass Tolerance:").pack()
tolerance_var = tk.StringVar(value="0.5")
tk.Entry(root, textvariable=tolerance_var).pack()

tk.Label(root, text="Stuctural Configurations").pack(pady=5)
branching_var = tk.BooleanVar(value=True)
tk.Checkbutton(root, text="Allow Ramifications", variable=branching_var).pack()
cycles_var = tk.BooleanVar(value=True)
tk.Checkbutton(root, text="Allow Cyclic", variable=cycles_var).pack()
double_var = tk.BooleanVar(value=True)
tk.Checkbutton(root, text="Allow Double Bonds", variable=double_var).pack()
triple_var = tk.BooleanVar(value=False)
tk.Checkbutton(root, text="Allow Triple Bonds", variable=triple_var).pack()

tk.Label(root, text="Maximum Double Bonds:").pack()
max_double_var = tk.StringVar(value="3")
tk.Entry(root, textvariable=max_double_var).pack()

tk.Label(root, text="Maxium Triple Bonds:").pack()
max_triple_var = tk.StringVar(value="2")
tk.Entry(root, textvariable=max_triple_var).pack()

tk.Label(root, text="Output File").pack(pady=5)
save_xyz_var = tk.BooleanVar(value=True)
tk.Checkbutton(root, text="Save .xyz", variable=save_xyz_var).pack()
save_mol_var = tk.BooleanVar(value=True)
tk.Checkbutton(root, text="Save .mol", variable=save_mol_var).pack()
save_svg_var = tk.BooleanVar(value=True)
tk.Checkbutton(root, text="Save .svg", variable=save_svg_var).pack()

tk.Button(root, text="Execute", command=run_program).pack(pady=10)
tk.Button(root, text="Open Output Folder", command=open_output_folder).pack()

root.mainloop()
