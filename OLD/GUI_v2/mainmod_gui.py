#################################################################################
###                                                                           ###
###                   ForMileS - FORMATION OF MASS SMILES                     ###
###                         (GUI Version - Cross-Platform)                    ###
###                                                                           ###
#################################################################################

import os
import sys
import re
from PyQt5.QtWidgets import (QApplication, QMainWindow, QTabWidget, QWidget, 
                            QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, 
                            QPushButton, QTextEdit, QFileDialog, QProgressBar,
                            QCheckBox, QSpinBox, QDoubleSpinBox, QListWidget,
                            QScrollArea, QGroupBox, QInputDialog, QMessageBox)
from PyQt5.QtCore import Qt, QThread, pyqtSignal
from PyQt5.QtGui import QPixmap
from submod1_gui import generate_smiles
from submod2_gui import filter_smiles
from submod3_gui import generate_charged_smiles, filter_charged_smiles_by_mass
from submod4_gui import smiles_to_molecules

class WorkerThread(QThread):
    progress_signal = pyqtSignal(int, str)
    finished_signal = pyqtSignal(bool)
    
    def __init__(self, function, *args, **kwargs):
        super().__init__()
        self.function = function
        self.args = args
        self.kwargs = kwargs
        self._is_running = True
    
    def run(self):
        try:
            self.function(*self.args, **self.kwargs)
            self.finished_signal.emit(True)
        except Exception as e:
            self.progress_signal.emit(100, f"Error: {str(e)}")
            self.finished_signal.emit(False)
        finally:
            self._is_running = False
    
    def stop(self):
        self._is_running = False
        self.wait(2000)  # Wait up to 2 seconds for clean exit
        if self.isRunning():
            self.terminate()  # Force termination if needed

class ForMileSGUI(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("ForMileS - Formation of Mass SMILES")
        self.setGeometry(100, 100, 1000, 800)
        
        self.central_widget = QTabWidget()
        self.setCentralWidget(self.central_widget)
        
        # Initialize variables
        self.output_dir = ""
        self.current_step = 0
        self.molecules = []
        self.data = []
        self.workers = []  # Track all worker threads
        
        # Create tabs
        self.setup_input_tab()
        self.setup_process_tab()
        self.setup_results_tab()
        
    def setup_input_tab(self):
        tab = QWidget()
        layout = QVBoxLayout()
        
        # Formula input
        formula_group = QGroupBox("Molecular Formula")
        formula_layout = QHBoxLayout()
        formula_layout.addWidget(QLabel("Formula:"))
        self.formula_input = QLineEdit("C6O2")
        formula_layout.addWidget(self.formula_input)
        formula_group.setLayout(formula_layout)
        
        # Charge input
        charge_group = QGroupBox("Charge")
        charge_layout = QHBoxLayout()
        charge_layout.addWidget(QLabel("Charge:"))
        self.charge_input = QSpinBox()
        self.charge_input.setRange(-5, 5)
        self.charge_input.setValue(1)
        charge_layout.addWidget(self.charge_input)
        charge_group.setLayout(charge_layout)
        
        # Precursor features
        precursor_group = QGroupBox("Precursor Features")
        precursor_layout = QVBoxLayout()
        self.precursor_list = QListWidget()
        self.precursor_list.addItem("CCCOCCCO")
        precursor_layout.addWidget(self.precursor_list)
        
        add_btn = QPushButton("Add Feature")
        add_btn.clicked.connect(self.add_precursor_feature)
        remove_btn = QPushButton("Remove Selected")
        remove_btn.clicked.connect(self.remove_precursor_feature)
        
        btn_layout = QHBoxLayout()
        btn_layout.addWidget(add_btn)
        btn_layout.addWidget(remove_btn)
        precursor_layout.addLayout(btn_layout)
        precursor_group.setLayout(precursor_layout)
        
        # Structure options
        structure_group = QGroupBox("Structure Options")
        structure_layout = QVBoxLayout()
        self.branched_check = QCheckBox("Allow Branched Structures")
        self.ring_check = QCheckBox("Allow Rings")
        structure_layout.addWidget(self.branched_check)
        structure_layout.addWidget(self.ring_check)
        structure_group.setLayout(structure_layout)
        
        # Mass filter
        mass_group = QGroupBox("Mass Filter")
        mass_layout = QVBoxLayout()
        
        target_layout = QHBoxLayout()
        target_layout.addWidget(QLabel("Target Mass:"))
        self.mass_input = QDoubleSpinBox()
        self.mass_input.setDecimals(4)
        self.mass_input.setRange(0, 10000)
        self.mass_input.setValue(117.092)
        target_layout.addWidget(self.mass_input)
        
        tolerance_layout = QHBoxLayout()
        tolerance_layout.addWidget(QLabel("Tolerance:"))
        self.tolerance_input = QDoubleSpinBox()
        self.tolerance_input.setDecimals(4)
        self.tolerance_input.setRange(0, 10)
        self.tolerance_input.setValue(0.5)
        tolerance_layout.addWidget(self.tolerance_input)
        
        mass_layout.addLayout(target_layout)
        mass_layout.addLayout(tolerance_layout)
        mass_group.setLayout(mass_layout)
        
        # Output directory
        output_group = QGroupBox("Output Directory")
        output_layout = QHBoxLayout()
        self.output_label = QLabel("Output will be saved in: [auto-generated]")
        browse_btn = QPushButton("Browse...")
        browse_btn.clicked.connect(self.select_output_dir)
        output_layout.addWidget(self.output_label)
        output_layout.addWidget(browse_btn)
        output_group.setLayout(output_layout)
        
        # Add all groups to layout
        layout.addWidget(formula_group)
        layout.addWidget(charge_group)
        layout.addWidget(precursor_group)
        layout.addWidget(structure_group)
        layout.addWidget(mass_group)
        layout.addWidget(output_group)
        layout.addStretch()
        
        tab.setLayout(layout)
        self.central_widget.addTab(tab, "Input Parameters")
    
    def add_precursor_feature(self):
        text, ok = QInputDialog.getText(self, "Add Precursor Feature", 
                                       "Enter SMILES/SMARTS pattern:")
        if ok and text:
            self.precursor_list.addItem(text)
    
    def remove_precursor_feature(self):
        for item in self.precursor_list.selectedItems():
            self.precursor_list.takeItem(self.precursor_list.row(item))
    
    def select_output_dir(self):
        dir_path = QFileDialog.getExistingDirectory(self, "Select Output Directory")
        if dir_path:
            self.output_dir = dir_path
            self.output_label.setText(f"Output will be saved in: {dir_path}")
    
    def setup_process_tab(self):
        tab = QWidget()
        layout = QVBoxLayout()
        
        # Progress bars
        self.progress_bars = []
        steps = [
            "1. Generate neutral SMILES",
            "2. Filter by precursor features",
            "3. Generate charged SMILES",
            "4. Filter by mass",
            "5. Generate molecule images"
        ]
        
        for step in steps:
            group = QGroupBox(step)
            step_layout = QVBoxLayout()
            progress = QProgressBar()
            status = QLabel("Ready")
            step_layout.addWidget(progress)
            step_layout.addWidget(status)
            group.setLayout(step_layout)
            layout.addWidget(group)
            self.progress_bars.append((progress, status))
        
        # Start button
        self.start_btn = QPushButton("Start Processing")
        self.start_btn.clicked.connect(self.start_processing)
        layout.addWidget(self.start_btn)
        
        # Console output
        self.console = QTextEdit()
        self.console.setReadOnly(True)
        layout.addWidget(self.console)
        
        tab.setLayout(layout)
        self.central_widget.addTab(tab, "Processing")
    
    def setup_results_tab(self):
        tab = QWidget()
        layout = QVBoxLayout()
        
        # Results display
        self.results_area = QScrollArea()
        self.results_container = QWidget()
        self.results_layout = QVBoxLayout()
        self.results_container.setLayout(self.results_layout)
        self.results_area.setWidget(self.results_container)
        self.results_area.setWidgetResizable(True)
        layout.addWidget(self.results_area)
        
        # Export button
        export_btn = QPushButton("Export Results")
        export_btn.clicked.connect(self.export_results)
        layout.addWidget(export_btn)
        
        tab.setLayout(layout)
        self.central_widget.addTab(tab, "Results")
    
    def log_message(self, message):
        self.console.append(message)
    
    def update_progress(self, step, value, message):
        progress, status = self.progress_bars[step]
        progress.setValue(value)
        status.setText(message)
        
        if value == 100:
            status.setStyleSheet("color: green;")
        elif value > 0:
            status.setStyleSheet("color: blue;")
        else:
            status.setStyleSheet("")
    
    def start_processing(self):
        # Get input values
        formula = self.formula_input.text().strip()
        charge = self.charge_input.value()
        precursor_features = [self.precursor_list.item(i).text() 
                            for i in range(self.precursor_list.count())]
        branched = self.branched_check.isChecked()
        ring = self.ring_check.isChecked()
        target_mass = self.mass_input.value()
        tolerance = self.tolerance_input.value()
        
        # Validate inputs
        if not formula:
            QMessageBox.critical(self, "Error", "Please enter a molecular formula")
            return
        
        if not precursor_features:
            reply = QMessageBox.question(self, "Warning", 
                                       "No precursor features specified. Continue anyway?",
                                       QMessageBox.Yes | QMessageBox.No)
            if reply == QMessageBox.No:
                return
        
        # Create output directory if not specified
        if not self.output_dir:
            self.output_dir = f"OutputFiles_{formula}_Charge_{charge}"
            os.makedirs(self.output_dir, exist_ok=True)
            self.output_label.setText(f"Output will be saved in: {os.path.abspath(self.output_dir)}")
        
        # Reset UI
        self.current_step = 0
        self.molecules = []
        self.data = []
        
        for progress, status in self.progress_bars:
            progress.setValue(0)
            status.setText("Ready")
            status.setStyleSheet("")
        
        self.console.clear()
        self.clear_results_tab()
        
        # Disable UI during processing
        self.setEnabled(False)
        self.start_btn.setText("Processing...")
        
        # Execute first step
        self.execute_step(0, generate_smiles, formula, charge, "nSMILES.txt")
    
    def execute_step(self, step, function, *args):
        worker = WorkerThread(function, *args)
        self.workers.append(worker)
        
        worker.progress_signal.connect(
            lambda p, m: self.update_progress(step, p, m))
        worker.finished_signal.connect(
            lambda success: self.on_step_finished(step, success))
        worker.finished_signal.connect(
            lambda: self.workers.remove(worker))
        
        worker.start()
    
    def on_step_finished(self, step, success):
        if not success:
            self.cleanup_after_processing()
            return
        
        # Get current parameters
        formula = self.formula_input.text().strip()
        charge = self.charge_input.value()
        precursor_features = [self.precursor_list.item(i).text() 
                            for i in range(self.precursor_list.count())]
        branched = self.branched_check.isChecked()
        ring = self.ring_check.isChecked()
        target_mass = self.mass_input.value()
        tolerance = self.tolerance_input.value()
        
        # Determine next step
        if step == 0:  # Neutral SMILES generated
            self.execute_step(1, filter_smiles, "nSMILES.txt", formula, charge, 
                            precursor_features, branched, ring, "ParentRelatedSMILES.txt")
        elif step == 1:  # Filtered by precursor
            self.execute_step(2, generate_charged_smiles, formula, charge, 
                            "ParentRelatedSMILES.txt", "chargedSMILES.txt")
        elif step == 2:  # Charged SMILES generated
            self.execute_step(3, filter_charged_smiles_by_mass, formula, charge, 
                            "chargedSMILES.txt", target_mass, tolerance, "filteredchargedSMILES.txt")
        elif step == 3:  # Filtered by mass
            self.execute_step(4, self.generate_and_display_molecules, formula, charge, 
                            "filteredchargedSMILES.txt")
        elif step == 4:  # All done
            self.log_message("\nProcessing completed successfully!")
            self.cleanup_after_processing(success=True)
    
    def generate_and_display_molecules(self, formula, charge, input_file):
        input_path = os.path.join(self.output_dir, input_file)
        
        try:
            self.molecules, self.data = smiles_to_molecules(formula, charge, input_path)
            self.display_results()
            self.update_progress(4, 100, f"Generated {len(self.molecules)} molecules")
            return True
        except Exception as e:
            self.log_message(f"\nError generating molecules: {str(e)}")
            self.update_progress(4, 100, f"Error: {str(e)}")
            return False
    
    def display_results(self):
        self.clear_results_tab()
        
        if not self.data:
            no_results = QLabel("No molecules were generated or an error occurred.")
            no_results.setAlignment(Qt.AlignCenter)
            self.results_layout.addWidget(no_results)
            return
        
        for idx, (smiles, mass, formula) in enumerate(self.data):
            group = QGroupBox(f"Molecule {idx + 1}")
            layout = QHBoxLayout()
            
            # Load and display image
            safe_smiles = self.sanitize_filename(smiles)
            img_path = os.path.join(self.output_dir, f"mol_{idx + 1}_{safe_smiles}.png")
            
            if os.path.exists(img_path):
                pixmap = QPixmap(img_path)
                img_label = QLabel()
                img_label.setPixmap(pixmap)
                img_label.setAlignment(Qt.AlignCenter)
                layout.addWidget(img_label)
            
            # Display info
            info_layout = QVBoxLayout()
            info_layout.addWidget(QLabel(f"<b>Formula:</b> {formula}"))
            info_layout.addWidget(QLabel(f"<b>Mass:</b> {mass}"))
            info_layout.addWidget(QLabel(f"<b>SMILES:</b> {smiles}"))
            info_layout.addStretch()
            layout.addLayout(info_layout)
            
            group.setLayout(layout)
            self.results_layout.addWidget(group)
        
        self.results_layout.addStretch()
    
    def clear_results_tab(self):
        # Clear previous results
        while self.results_layout.count():
            item = self.results_layout.takeAt(0)
            widget = item.widget()
            if widget:
                widget.deleteLater()
    
    def export_results(self):
        if not self.data:
            QMessageBox.information(self, "Export", "No results to export")
            return
        
        file_path, _ = QFileDialog.getSaveFileName(
            self, "Export Results", "", "CSV Files (*.csv);;All Files (*)")
        
        if file_path:
            try:
                with open(file_path, 'w') as f:
                    f.write("Index,Formula,Mass,SMILES\n")
                    for idx, (smiles, mass, formula) in enumerate(self.data):
                        f.write(f"{idx+1},{formula},{mass},{smiles}\n")
                QMessageBox.information(self, "Export", "Results exported successfully")
            except Exception as e:
                QMessageBox.critical(self, "Export Error", f"Failed to export: {str(e)}")
    
    def cleanup_after_processing(self, success=False):
        # Re-enable UI
        self.setEnabled(True)
        self.start_btn.setText("Start Processing")
        
        if success:
            self.central_widget.setCurrentIndex(2)  # Switch to results tab
    
    def closeEvent(self, event):
        """Clean up threads when window is closed"""
        running_threads = False
        
        for worker in self.workers:
            if worker.isRunning():
                running_threads = True
                worker.stop()
        
        if running_threads:
            reply = QMessageBox.question(
                self, "Threads Running",
                "Background processes are still running. Force quit?",
                QMessageBox.Yes | QMessageBox.No)
            
            if reply == QMessageBox.No:
                event.ignore()
                return
        
        event.accept()
    
    def sanitize_filename(self, smiles):
        return re.sub(r'[^a-zA-Z0-9._-]', '_', smiles)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    
    # Set some application-wide styles
    app.setStyle("Fusion")
    
    window = ForMileSGUI()
    window.show()
    
    try:
        sys.exit(app.exec_())
    except Exception as e:
        print(f"Application error: {str(e)}")