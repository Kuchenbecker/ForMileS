"""
main.py - Entry point for ForMileS GUI app.
Use with PyInstaller to create a standalone executable.

Build with:
    pyinstaller --onefile --windowed main.py
"""

import sys
import os

# If needed, adjust paths when bundled
def get_resource_path(relative_path):
    """ Get absolute path to resource, works for dev and PyInstaller """
    try:
        # PyInstaller creates a temp folder and stores path in _MEIPASS
        base_path = sys._MEIPASS
    except AttributeError:
        base_path = os.path.abspath(".")

    return os.path.join(base_path, relative_path)

# Optional: Fix matplotlib backend if needed by RDKit
import matplotlib
matplotlib.use('Agg')  # use non-GUI backend

# Start the main GUI application
import gui_formiles

# gui_formiles.py includes `root.mainloop()` so no further call is needed here
