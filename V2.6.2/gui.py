import tkinter as tk
from tkinter import ttk, messagebox
import json
import subprocess
import os
import sys
from threading import Thread

CONFIG_FILE = "config.json"

# Add this function to handle resource paths
def resource_path(relative_path):
    """ Get absolute path to resource, works for dev and for PyInstaller """
    try:
        # PyInstaller creates a temp folder and stores path in _MEIPASS
        base_path = sys._MEIPASS
    except Exception:
        base_path = os.path.abspath(".")

    return os.path.join(base_path, relative_path)

def get_output_folder():
    try:
        with open(CONFIG_FILE) as f:
            config = json.load(f)
        formula = config.get("FORMULA", "Output")
        charge = config.get("CHARGE", 0)
        return f"OutputFiles_{formula}_Charge_{charge}"
    except:
        return "."

def show_waiting_window():
    """Create and display a waiting window"""
    waiting_window = tk.Toplevel()
    waiting_window.title("Processing")
    waiting_window.geometry("300x100")
    
    # Center the window
    waiting_window.update_idletasks()
    width = waiting_window.winfo_width()
    height = waiting_window.winfo_height()
    x = (waiting_window.winfo_screenwidth() // 2) - (width // 2)
    y = (waiting_window.winfo_screenheight() // 2) - (height // 2)
    waiting_window.geometry(f"+{x}+{y}")
    
    tk.Label(waiting_window, text="ForMileS is Generating Structures.\nPlease, wait...", 
             font=('Helvetica', 12)).pack(pady=20)
    
    # Make the window modal
    waiting_window.grab_set()
    waiting_window.transient(root)
    waiting_window.focus_force()
    
    # Prevent closing while processing
    waiting_window.protocol("WM_DELETE_WINDOW", lambda: None)
    
    return waiting_window

def run_program_thread():
    """Thread function to run the program and manage windows"""
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
        
        # Close waiting window and show success message
        root.after(0, lambda: [waiting_window.destroy(), 
                              messagebox.showinfo("Finished", "ForMileS Ended Successfully. Check output directory.")])
    except Exception as e:
        # Close waiting window and show error message
        root.after(0, lambda: [waiting_window.destroy(), 
                              messagebox.showerror("Error", f"An error occurred trying to execute: {e}")])

def run_program():
    """Start the program in a separate thread and show waiting window"""
    global waiting_window
    waiting_window = show_waiting_window()
    # Start the program in a separate thread to avoid freezing the GUI
    Thread(target=run_program_thread, daemon=True).start()

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
root.title("ForMileS V2.6")
root.geometry("400x630")

tk.Label(root, text="Molecular Formula:").pack()
formula_var = tk.StringVar(value="C6O2")
tk.Entry(root, textvariable=formula_var).pack()

tk.Label(root, text="Basic Molecular Scaffold (BMS):").pack()
scaffold_var = tk.StringVar(value="CCCOC")
tk.Entry(root, textvariable=scaffold_var).pack()

tk.Label(root, text="Charge:").pack()
charge_var = tk.StringVar(value="+1")
tk.Entry(root, textvariable=charge_var).pack()

tk.Label(root, text="Exact Mass:").pack()
mass_var = tk.StringVar(value="117.092")
tk.Entry(root, textvariable=mass_var).pack()

tk.Label(root, text="Mass Tolerance:").pack()
tolerance_var = tk.StringVar(value="0.5")
tk.Entry(root, textvariable=tolerance_var).pack()

tk.Label(root, text="Stuctural Configurations").pack(pady=5)
branching_var = tk.BooleanVar(value=False)
tk.Checkbutton(root, text="Allow Ramifications", variable=branching_var).pack()
cycles_var = tk.BooleanVar(value=False)
tk.Checkbutton(root, text="Allow Cyclic", variable=cycles_var).pack()
double_var = tk.BooleanVar(value=False)
tk.Checkbutton(root, text="Allow Double Bonds", variable=double_var).pack()
triple_var = tk.BooleanVar(value=False)
tk.Checkbutton(root, text="Allow Triple Bonds", variable=triple_var).pack()

tk.Label(root, text="Maximum Double Bonds:").pack()
max_double_var = tk.StringVar(value="1")
tk.Entry(root, textvariable=max_double_var).pack()

tk.Label(root, text="Maxium Triple Bonds:").pack()
max_triple_var = tk.StringVar(value="1")
tk.Entry(root, textvariable=max_triple_var).pack()

tk.Label(root, text="Output File").pack(pady=5)
save_xyz_var = tk.BooleanVar(value=False)
tk.Checkbutton(root, text="Save .xyz", variable=save_xyz_var).pack()
save_mol_var = tk.BooleanVar(value=False)
tk.Checkbutton(root, text="Save .mol", variable=save_mol_var).pack()
save_svg_var = tk.BooleanVar(value=True)
tk.Checkbutton(root, text="Save .svg", variable=save_svg_var).pack()

tk.Button(root, text="Execute", command=run_program).pack(pady=10)
tk.Button(root, text="Open Output Folder", command=open_output_folder).pack()

root.mainloop()