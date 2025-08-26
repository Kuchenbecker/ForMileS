<<<<<<< HEAD
# Updated gui_formiles.py with title, description, logo, runtime, kill button, and close confirmation
=======
>>>>>>> 965c5a91956b851a1a46d64358e56db8b3aae644
import tkinter as tk
from tkinter import ttk, messagebox
import threading
import os
<<<<<<< HEAD
import time
from PIL import ImageTk, Image
import sys
=======
import multiprocessing
from PIL import ImageTk, Image
>>>>>>> 965c5a91956b851a1a46d64358e56db8b3aae644

from submod1 import generate_smiles
from submod2 import filter_smiles
from submod3 import generate_charged_smiles, filter_charged_smiles_by_mass
from submod4 import smiles_to_molecules

<<<<<<< HEAD
# Global variables
generated_images = []
formiles_thread = None
stop_requested = False

# GUI setup
root = tk.Tk()
root.title("ForMileS - Formation of Mass SMILES")

main_frame = ttk.Frame(root)
main_frame.pack(fill="both", expand=True)

left_frame = ttk.Frame(main_frame)
left_frame.pack(side="left", padx=10, pady=10)

right_frame = ttk.Frame(main_frame)
right_frame.pack(side="right", fill="both", expand=True)

img_canvas = tk.Canvas(right_frame)
img_scrollbar = ttk.Scrollbar(right_frame, orient="vertical", command=img_canvas.yview)
img_canvas.configure(yscrollcommand=img_scrollbar.set)
img_scrollbar.pack(side="right", fill="y")
img_canvas.pack(side="left", fill="both", expand=True)
img_frame = ttk.Frame(img_canvas)
img_canvas.create_window((0, 0), window=img_frame, anchor="nw")

progress_var = tk.DoubleVar()
progress_label = ttk.Label(left_frame, text="Progress: 0%")
status_label = ttk.Label(left_frame, text="Status: Idle")
runtime_label = ttk.Label(left_frame, text="Runtime: 0.0 s")

# Load logo
try:
    logo_img = Image.open("logo.png")
    logo_img = logo_img.resize((50, 50))
    logo_photo = ImageTk.PhotoImage(logo_img)
    logo_label = ttk.Label(left_frame, image=logo_photo)
    logo_label.image = logo_photo
    logo_label.grid(row=0, column=0, rowspan=2, padx=(0, 10), pady=(5, 5))
except Exception as e:
    print(f"Logo not loaded: {e}")

# Title and description
title_label = ttk.Label(left_frame, text="ForMileS", font=("Helvetica", 20, "bold"))
title_label.grid(row=0, column=1, sticky="w", pady=(5, 0))
desc_label = ttk.Label(left_frame, text="Formation of Mass SMILES GUI Tool", font=("Helvetica", 10))
desc_label.grid(row=1, column=1, sticky="w", pady=(0, 15))

def update_progress(value, status_text=""):
    progress_var.set(value)
    progress_label.config(text=f"Progress: {int(value)}%")
    if status_text:
        status_label.config(text=f"Status: {status_text}")
    root.update_idletasks()

def run_formiles(as_svg=False):
    global stop_requested
    stop_requested = False
    start_time = time.time()
=======
# Global variable to store generated image paths
generated_images = []

def update_progress(value):
    progress_var.set(value)
    root.update_idletasks()

def run_formiles(as_svg=False):
>>>>>>> 965c5a91956b851a1a46d64358e56db8b3aae644
    try:
        formula = formula_entry.get().strip()
        charge = int(charge_entry.get().strip())
        target_mass = float(mass_entry.get().strip())
        features = features_text.get("1.0", "end").strip().split("\n")
        is_branched = branched_var.get()
        has_ring = ring_var.get()
<<<<<<< HEAD
=======
        try:
            num_processes = min(max(1, int(cores_var.get())), multiprocessing.cpu_count())
        except:
            num_processes = 1
>>>>>>> 965c5a91956b851a1a46d64358e56db8b3aae644

        output_dir = f"OutputFiles_{formula}_Charge_{charge}"
        os.makedirs(output_dir, exist_ok=True)

<<<<<<< HEAD
        update_progress(5, "Generating SMILES")
        if stop_requested: return
        generate_smiles(formula, charge, f"nSMILES_{formula}.txt", should_stop=lambda: stop_requested, update_callback=gui_progress_callback)

        update_progress(25, "Filtering SMILES")
        if stop_requested: return
        filter_smiles(f"nSMILES_{formula}.txt", formula, charge, features, is_branched, has_ring, f"ParentRelatedSMILES_{formula}.txt", should_stop=lambda: stop_requested, update_callback=gui_progress_callback)

        update_progress(50, "Generating charged SMILES")
        if stop_requested: return
        generate_charged_smiles(formula, charge, f"ParentRelatedSMILES_{formula}.txt", f"chargedSMILES_{formula}.txt", should_stop=lambda: stop_requested, update_callback=gui_progress_callback)

        update_progress(70, "Filtering charged by mass")
        if stop_requested: return
        filter_charged_smiles_by_mass(formula, charge, f"chargedSMILES_{formula}.txt", target_mass, 0.05, f"filteredchargedSMILES_{formula}.txt", should_stop=lambda: stop_requested, update_callback=gui_progress_callback)

        update_progress(90, "Generating images")
        if stop_requested: return
        global generated_images
        _, _, generated_images = smiles_to_molecules(formula, charge, f"filteredchargedSMILES_{formula}.txt", as_svg=as_svg)

        update_progress(100, "Finished")
        elapsed = time.time() - start_time
        runtime_label.config(text=f"Runtime: {elapsed:.2f} s")
        messagebox.showinfo("Done", "ForMileS completed successfully.\nCheck the OutputFiles folder.")
    except Exception as e:
        messagebox.showerror("Error", f"Something went wrong:\n{str(e)}")
        update_progress(0, "Idle")
        runtime_label.config(text="Runtime: --")

def gui_progress_callback(current, total):
    percent = int((current / total) * 20)  # adjust scaling for total progress
    update_progress(percent)

def run_formiles_threaded():
    global formiles_thread
    formiles_thread = threading.Thread(target=run_formiles, args=(svg_var.get(),))
    formiles_thread.start()



def kill_processing():
    global stop_requested
    stop_requested = True
    status_label.config(text="Status: Aborting...")

# Ask confirmation on close if thread is running
def on_closing():
    if formiles_thread and formiles_thread.is_alive():
        if messagebox.askyesno("Confirm Exit", "ForMileS is still running. Are you sure you want to close and stop processing?"):
            kill_processing()
            root.destroy()
    else:
        root.destroy()

root.protocol("WM_DELETE_WINDOW", on_closing)
=======
        update_progress(5)
        generate_smiles(formula, charge, f"nSMILES_{formula}.txt", num_processes)
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
>>>>>>> 965c5a91956b851a1a46d64358e56db8b3aae644

def show_images():
    if not generated_images:
        messagebox.showwarning("No Images", "Run ForMileS first to generate images.")
        return

<<<<<<< HEAD
    for widget in img_frame.winfo_children():
        widget.destroy()
=======
    viewer = tk.Toplevel(root)
    viewer.title("Generated Molecule Images")

    canvas = tk.Canvas(viewer)
    scrollbar = ttk.Scrollbar(viewer, orient="vertical", command=canvas.yview)
    canvas.configure(yscrollcommand=scrollbar.set)

    scrollbar.pack(side="right", fill="y")
    canvas.pack(side="left", fill="both", expand=True)

    frame = ttk.Frame(canvas)
    canvas.create_window((0, 0), window=frame, anchor="nw")
>>>>>>> 965c5a91956b851a1a46d64358e56db8b3aae644

    for img_path in generated_images:
        if img_path.endswith(".png"):
            try:
                img = Image.open(img_path)
                img.thumbnail((300, 300))
                photo = ImageTk.PhotoImage(img)
<<<<<<< HEAD
                label = ttk.Label(img_frame, image=photo)
                label.image = photo
=======
                label = ttk.Label(frame, image=photo)
                label.image = photo  # keep a reference!
>>>>>>> 965c5a91956b851a1a46d64358e56db8b3aae644
                label.pack(padx=5, pady=5)
            except Exception as e:
                print(f"Could not open image: {img_path}\n{e}")

<<<<<<< HEAD
    img_frame.update_idletasks()
    img_canvas.config(scrollregion=img_canvas.bbox("all"))

# Inputs start at row 2 because of title block
ttk.Label(left_frame, text="Molecular Formula (e.g., C4O2):").grid(row=2, column=0, columnspan=2, sticky="w")
formula_entry = ttk.Entry(left_frame, width=30)
formula_entry.grid(row=3, column=0, columnspan=2)

ttk.Label(left_frame, text="Charge (+1 or -1):").grid(row=4, column=0, columnspan=2, sticky="w")
charge_entry = ttk.Entry(left_frame, width=30)
charge_entry.grid(row=5, column=0, columnspan=2)

ttk.Label(left_frame, text="Target Mass (e.g., 89.060):").grid(row=6, column=0, columnspan=2, sticky="w")
mass_entry = ttk.Entry(left_frame, width=30)
mass_entry.grid(row=7, column=0, columnspan=2)

ttk.Label(left_frame, text="Precursor SMARTS (one per line):").grid(row=8, column=0, columnspan=2, sticky="w")
features_text = tk.Text(left_frame, width=30, height=5)
features_text.grid(row=9, column=0, columnspan=2)

=======
    frame.update_idletasks()
    canvas.config(scrollregion=canvas.bbox("all"))

# GUI setup
root = tk.Tk()
root.title("ForMileS - Formation of Mass SMILES")

# Create main frame
mainframe = ttk.Frame(root, padding="10")
mainframe.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))

# Input fields
ttk.Label(mainframe, text="Molecular Formula (e.g., C4O2):").grid(row=0, column=0, sticky="w", pady=5)
formula_entry = ttk.Entry(mainframe, width=30)
formula_entry.grid(row=0, column=1, pady=5)

ttk.Label(mainframe, text="Charge (+1 or -1):").grid(row=1, column=0, sticky="w", pady=5)
charge_entry = ttk.Entry(mainframe, width=30)
charge_entry.grid(row=1, column=1, pady=5)

ttk.Label(mainframe, text="Target Mass (e.g., 89.060):").grid(row=2, column=0, sticky="w", pady=5)
mass_entry = ttk.Entry(mainframe, width=30)
mass_entry.grid(row=2, column=1, pady=5)

ttk.Label(mainframe, text="Precursor SMARTS (one per line):").grid(row=3, column=0, sticky="nw", pady=5)
features_text = tk.Text(mainframe, width=30, height=5)
features_text.grid(row=3, column=1, pady=5)

# Checkboxes
>>>>>>> 965c5a91956b851a1a46d64358e56db8b3aae644
branched_var = tk.BooleanVar(value=False)
ring_var = tk.BooleanVar(value=False)
svg_var = tk.BooleanVar(value=False)

<<<<<<< HEAD
ttk.Checkbutton(left_frame, text="Allow Branched", variable=branched_var).grid(row=10, column=0, sticky="w")
ttk.Checkbutton(left_frame, text="Allow Rings", variable=ring_var).grid(row=10, column=1, sticky="w")
ttk.Checkbutton(left_frame, text="Generate SVG (PNG always generated)", variable=svg_var).grid(row=11, column=0, columnspan=2, sticky="w")

ttk.Button(left_frame, text="Run ForMileS", command=run_formiles_threaded).grid(row=12, column=0, columnspan=2, pady=10)
ttk.Button(left_frame, text="Kill Processing", command=kill_processing).grid(row=13, column=0, columnspan=2, pady=(0, 10))

ttk.Progressbar(left_frame, variable=progress_var, maximum=100).grid(row=14, column=0, columnspan=2, pady=10, sticky="ew")
progress_label.grid(row=15, column=0, columnspan=2)
status_label.grid(row=16, column=0, columnspan=2)
runtime_label.grid(row=17, column=0, columnspan=2)

ttk.Button(left_frame, text="Show Images Output", command=show_images).grid(row=18, column=0, columnspan=2, pady=5)

root.mainloop()
=======
ttk.Checkbutton(mainframe, text="Allow Branched", variable=branched_var).grid(row=4, column=0, sticky="w", pady=5)
ttk.Checkbutton(mainframe, text="Allow Rings", variable=ring_var).grid(row=4, column=1, sticky="w", pady=5)

# CPU cores selection
max_cores = multiprocessing.cpu_count()
ttk.Label(mainframe, text=f"CPU Cores (1-{max_cores}):").grid(row=5, column=0, sticky="w", pady=5)
cores_var = tk.StringVar(value=str(max(1, max_cores//2)))
cores_entry = ttk.Entry(mainframe, width=5, textvariable=cores_var)
cores_entry.grid(row=5, column=1, sticky="w", pady=5)

ttk.Checkbutton(mainframe, text="Generate SVG instead of PNG", variable=svg_var).grid(row=6, column=0, columnspan=2, sticky="w", pady=5)

# Buttons
ttk.Button(mainframe, text="Run ForMileS", command=run_formiles_threaded).grid(row=7, column=0, columnspan=2, pady=10)

# Progress bar
progress_var = tk.DoubleVar()
progress = ttk.Progressbar(mainframe, variable=progress_var, maximum=100)
progress.grid(row=8, column=0, columnspan=2, pady=10, sticky="ew")

# Show Images button
ttk.Button(mainframe, text="Show Images Output", command=show_images).grid(row=9, column=0, columnspan=2, pady=5)

root.mainloop()
>>>>>>> 965c5a91956b851a1a46d64358e56db8b3aae644
