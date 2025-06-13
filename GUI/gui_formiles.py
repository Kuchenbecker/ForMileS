import tkinter as tk
from tkinter import ttk, messagebox
import threading
import os
from PIL import ImageTk, Image

from submod1 import generate_smiles
from submod2 import filter_smiles
from submod3 import generate_charged_smiles, filter_charged_smiles_by_mass
from submod4 import smiles_to_molecules

# Global variable to store generated image paths
generated_images = []

def update_progress(value):
    progress_var.set(value)
    root.update_idletasks()

def run_formiles(as_svg=False):
    try:
        formula = formula_entry.get().strip()
        charge = int(charge_entry.get().strip())
        target_mass = float(mass_entry.get().strip())
        features = features_text.get("1.0", "end").strip().split("\n")
        is_branched = branched_var.get()
        has_ring = ring_var.get()

        output_dir = f"OutputFiles_{formula}_Charge_{charge}"
        os.makedirs(output_dir, exist_ok=True)

        update_progress(5)
        generate_smiles(formula, charge, f"nSMILES_{formula}.txt")
        update_progress(25)

        filter_smiles(f"nSMILES_{formula}.txt", formula, charge, features, is_branched, has_ring, f"ParentRelatedSMILES_{formula}.txt")
        update_progress(50)

        generate_charged_smiles(formula, charge, f"ParentRelatedSMILES_{formula}.txt", f"chargedSMILES_{formula}.txt")
        update_progress(70)

        filter_charged_smiles_by_mass(formula, charge, f"chargedSMILES_{formula}.txt", target_mass, 0.05, f"filteredchargedSMILES_{formula}.txt")
        update_progress(90)

        global generated_images
        _, _, generated_images = smiles_to_molecules(formula, charge, f"filteredchargedSMILES_{formula}.txt", as_svg=as_svg)

        update_progress(100)
        messagebox.showinfo("Done", "ForMileS completed successfully.\nCheck the OutputFiles folder.")
    except Exception as e:
        messagebox.showerror("Error", f"Something went wrong:\n{str(e)}")
        update_progress(0)

def run_formiles_threaded():
    thread = threading.Thread(target=run_formiles, args=(svg_var.get(),))
    thread.start()

def show_images():
    if not generated_images:
        messagebox.showwarning("No Images", "Run ForMileS first to generate images.")
        return

    viewer = tk.Toplevel(root)
    viewer.title("Generated Molecule Images")

    canvas = tk.Canvas(viewer)
    scrollbar = ttk.Scrollbar(viewer, orient="vertical", command=canvas.yview)
    canvas.configure(yscrollcommand=scrollbar.set)

    scrollbar.pack(side="right", fill="y")
    canvas.pack(side="left", fill="both", expand=True)

    frame = ttk.Frame(canvas)
    canvas.create_window((0, 0), window=frame, anchor="nw")

    for img_path in generated_images:
        if img_path.endswith(".png"):
            try:
                img = Image.open(img_path)
                img.thumbnail((300, 300))
                photo = ImageTk.PhotoImage(img)
                label = ttk.Label(frame, image=photo)
                label.image = photo  # keep a reference!
                label.pack(padx=5, pady=5)
            except Exception as e:
                print(f"Could not open image: {img_path}\n{e}")

    frame.update_idletasks()
    canvas.config(scrollregion=canvas.bbox("all"))

# GUI setup
root = tk.Tk()
root.title("ForMileS - Formation of Mass SMILES")

ttk.Label(root, text="Molecular Formula (e.g., C4O2):").grid(row=0, column=0, sticky="w")
formula_entry = ttk.Entry(root, width=30)
formula_entry.grid(row=0, column=1)

ttk.Label(root, text="Charge (+1 or -1):").grid(row=1, column=0, sticky="w")
charge_entry = ttk.Entry(root, width=30)
charge_entry.grid(row=1, column=1)

ttk.Label(root, text="Target Mass (e.g., 89.060):").grid(row=2, column=0, sticky="w")
mass_entry = ttk.Entry(root, width=30)
mass_entry.grid(row=2, column=1)

ttk.Label(root, text="Precursor SMARTS (one per line):").grid(row=3, column=0, sticky="nw")
features_text = tk.Text(root, width=30, height=5)
features_text.grid(row=3, column=1)

branched_var = tk.BooleanVar(value=False)
ring_var = tk.BooleanVar(value=False)
svg_var = tk.BooleanVar(value=False)

ttk.Checkbutton(root, text="Allow Branched", variable=branched_var).grid(row=4, column=0, sticky="w")
ttk.Checkbutton(root, text="Allow Rings", variable=ring_var).grid(row=4, column=1, sticky="w")
ttk.Checkbutton(root, text="Generate SVG instead of PNG", variable=svg_var).grid(row=5, column=0, columnspan=2, sticky="w")

ttk.Button(root, text="Run ForMileS", command=run_formiles_threaded).grid(row=6, column=0, columnspan=2, pady=10)

# Progress bar
progress_var = tk.DoubleVar()
progress = ttk.Progressbar(root, variable=progress_var, maximum=100)
progress.grid(row=7, column=0, columnspan=2, pady=10, sticky="ew")

# Show Images button
ttk.Button(root, text="Show Images Output", command=show_images).grid(row=8, column=0, columnspan=2, pady=5)

root.mainloop()
