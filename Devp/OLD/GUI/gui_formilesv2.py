import tkinter as tk
from tkinter import ttk, messagebox
import threading
import os
import time
from PIL import ImageTk, Image
import sys

from gui_submod1 import generate_smiles
from gui_submod2 import filter_smiles
from gui_submod3 import generate_charged_smiles, filter_charged_smiles_by_mass
from gui_submod4 import smiles_to_molecules

# Global variables
generated_images = []
formiles_thread = None
stop_requested = False
start_time = 0

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
eta_label = ttk.Label(left_frame, text="ETA: --")
rate_label = ttk.Label(left_frame, text="Speed: -- it/s")

def format_time(seconds):
    if seconds < 60:
        return f"{seconds:.0f}s"
    elif seconds < 3600:
        return f"{seconds//60:.0f}m {seconds%60:.0f}s"
    else:
        return f"{seconds//3600:.0f}h {(seconds%3600)//60:.0f}m"

def update_progress(value, status_text="", current=None, total=None):
    # Only update if significant change or new status
    if abs(value - progress_var.get()) > 1 or status_text != status_label.cget("text"):
        progress_var.set(value)
        progress_label.config(text=f"Progress: {int(value)}%")
        status_label.config(text=f"Status: {status_text}")
    
    # Update time metrics less frequently
    if current is None or current % 100 == 0:
        elapsed = time.time() - start_time
        runtime_label.config(text=f"Runtime: {format_time(elapsed)}")
        
        if current is not None and total is not None and current > 0:
            rate = current / elapsed if elapsed > 0 else 0
            remaining = (total - current) / rate if rate > 0 else 0
            eta_label.config(text=f"ETA: {format_time(remaining)}")
            rate_label.config(text=f"Speed: {rate:.1f} it/s")
    
    # Throttle UI updates
    if current is None or current % 10 == 0:
        root.update_idletasks()

def run_formiles(as_svg=False):
    global stop_requested, start_time
    stop_requested = False
    start_time = time.time()
    
    # Disable GUI updates during heavy computation
    root.config(cursor="watch")
    root.update()
    
    eta_label.config(text="ETA: --")
    rate_label.config(text="Speed: -- it/s")
    
    try:
        formula = formula_entry.get().strip()
        charge = int(charge_entry.get().strip())
        target_mass = float(mass_entry.get().strip())
        features = features_text.get("1.0", "end").strip().split("\n")
        is_branched = branched_var.get()
        has_ring = ring_var.get()

        output_dir = f"OutputFiles_{formula}_Charge_{charge}"
        os.makedirs(output_dir, exist_ok=True)

        update_progress(5, "Generating SMILES")
        if stop_requested: return
        generate_smiles(formula, charge, f"nSMILES_{formula}.txt", 
                       should_stop=lambda: stop_requested, 
                       update_callback=lambda c, t: update_progress(
                           5 + 20*c/t, "Generating SMILES", c, t))

        update_progress(25, "Filtering SMILES")
        if stop_requested: return
        filter_smiles(f"nSMILES_{formula}.txt", formula, charge, features, 
                     is_branched, has_ring, f"ParentRelatedSMILES_{formula}.txt", 
                     should_stop=lambda: stop_requested,
                     update_callback=lambda c, t: update_progress(
                         25 + 25*c/t, "Filtering SMILES", c, t))

        update_progress(50, "Generating charged SMILES")
        if stop_requested: return
        generate_charged_smiles(formula, charge, f"ParentRelatedSMILES_{formula}.txt", 
                              f"chargedSMILES_{formula}.txt", 
                              should_stop=lambda: stop_requested,
                              update_callback=lambda c, t: update_progress(
                                  50 + 20*c/t, "Generating charged SMILES", c, t))

        update_progress(70, "Filtering by mass")
        if stop_requested: return
        filter_charged_smiles_by_mass(formula, charge, f"chargedSMILES_{formula}.txt", 
                                    target_mass, 0.05, f"filteredchargedSMILES_{formula}.txt", 
                                    should_stop=lambda: stop_requested,
                                    update_callback=lambda c, t: update_progress(
                                        70 + 20*c/t, "Filtering by mass", c, t))

        update_progress(90, "Generating images")
        if stop_requested: return
        global generated_images
        _, _, generated_images = smiles_to_molecules(formula, charge, 
                                                  f"filteredchargedSMILES_{formula}.txt", 
                                                  as_svg=as_svg)

        update_progress(100, "Finished")
        messagebox.showinfo("Done", "ForMileS completed successfully.\nCheck the OutputFiles folder.")
    except Exception as e:
        messagebox.showerror("Error", f"Something went wrong:\n{str(e)}")
        update_progress(0, "Idle")
    finally:
        # Re-enable GUI
        root.config(cursor="")
        root.update()

def run_formiles_threaded():
    global formiles_thread, stop_requested
    stop_requested = False
    formiles_thread = threading.Thread(target=run_formiles, args=(svg_var.get(),))
    formiles_thread.start()

def kill_processing():
    global stop_requested
    stop_requested = True
    status_label.config(text="Status: Aborting...")
    if formiles_thread and formiles_thread.is_alive():
        formiles_thread.join(timeout=1.0)

def on_closing():
    kill_processing()
    root.destroy()

root.protocol("WM_DELETE_WINDOW", on_closing)

def show_images():
    if not generated_images:
        messagebox.showwarning("No Images", "Run ForMileS first to generate images.")
        return

    for widget in img_frame.winfo_children():
        widget.destroy()

    for img_path in generated_images:
        if img_path.endswith(".png"):
            try:
                img = Image.open(img_path)
                img.thumbnail((300, 300))
                photo = ImageTk.PhotoImage(img)
                label = ttk.Label(img_frame, image=photo)
                label.image = photo
                label.pack(padx=5, pady=5)
            except Exception as e:
                print(f"Could not open image: {img_path}\n{e}")

    img_frame.update_idletasks()
    img_canvas.config(scrollregion=img_canvas.bbox("all"))

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

# Input fields
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

# Checkboxes
branched_var = tk.BooleanVar(value=False)
ring_var = tk.BooleanVar(value=False)
svg_var = tk.BooleanVar(value=False)

ttk.Checkbutton(left_frame, text="Allow Branched", variable=branched_var).grid(row=10, column=0, sticky="w")
ttk.Checkbutton(left_frame, text="Allow Rings", variable=ring_var).grid(row=10, column=1, sticky="w")
ttk.Checkbutton(left_frame, text="Generate SVG (PNG always generated)", variable=svg_var).grid(row=11, column=0, columnspan=2, sticky="w")

# Buttons
ttk.Button(left_frame, text="Run ForMileS", command=run_formiles_threaded).grid(row=12, column=0, columnspan=2, pady=10)
ttk.Button(left_frame, text="Kill Processing", command=kill_processing).grid(row=13, column=0, columnspan=2, pady=(0, 10))

# Progress area
ttk.Progressbar(left_frame, variable=progress_var, maximum=100).grid(row=14, column=0, columnspan=2, pady=10, sticky="ew")
progress_label.grid(row=15, column=0, columnspan=2)
status_label.grid(row=16, column=0, columnspan=2)
runtime_label.grid(row=17, column=0, columnspan=2)
eta_label.grid(row=18, column=0, columnspan=2)
rate_label.grid(row=19, column=0, columnspan=2)

ttk.Button(left_frame, text="Show Images Output", command=show_images).grid(row=20, column=0, columnspan=2, pady=5)

root.mainloop()