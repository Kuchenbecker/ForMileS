#################################################################################
###                                                                           ###
###                   ForMileS - FORMATION OF MASS SMILES                     ###
###                         (GUI Version - Cross-Platform)                    ###
###                                                                           ###
#################################################################################

import os
import sys
import re
import math
from PyQt5.QtWidgets import (QApplication, QMainWindow, QTabWidget, QWidget, 
                            QVBoxLayout, QHBoxLayout, QLabel, QLineEdit, 
                            QPushButton, QTextEdit, QFileDialog, QProgressBar,
                            QCheckBox, QSpinBox, QDoubleSpinBox, QListWidget,
                            QScrollArea, QGroupBox, QInputDialog, QMessageBox)
from PyQt5.QtCore import Qt, QThread, pyqtSignal
from PyQt5.QtGui import QPixmap
import networkx as nx
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors, rdMolDescriptors
from itertools import combinations, product
from PIL import Image, ImageDraw, ImageFont

class WorkerThread(QThread):
    progress_signal = pyqtSignal(int, str)
    finished_signal = pyqtSignal(bool)
    
    def __init__(self, function, *args, **kwargs):
        super().__init__()
        self.function = function
        self.args = args
        self.kwargs = kwargs
        self._is_running = True
        self.kwargs['progress_callback'] = self.progress_signal.emit
    
    def run(self):
        try:
            result = self.function(*self.args, **self.kwargs)
            self.finished_signal.emit(True)
        except Exception as e:
            self.progress_signal.emit(100, f"Error: {str(e)}")
            self.finished_signal.emit(False)
        finally:
            self._is_running = False
    
    def stop(self):
        self._is_running = False
        self.wait(2000)
        if self.isRunning():
            self.terminate()

class ForMileSGUI(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("ForMileS - Formation of Mass SMILES")
        self.setGeometry(100, 100, 1000, 800)
        
        self.central_widget = QTabWidget()
        self.setCentralWidget(self.central_widget)
        
        self.output_dir = ""
        self.current_step = 0
        self.molecules = []
        self.data = []
        self.workers = []
        self.current_step_results = {}
        
        # Chemistry configuration
        self.bond_orders = {
            ("C", "C"): [1, 2],
            ("C", "O"): [1, 2],
            ("O", "O"): [1],
        }
        self.max_valence = {
            "C": 4,
            "O": 2,
        }
        
        self.setup_ui()
    
    def setup_ui(self):
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
        self.current_step_results = {}
        
        for progress, status in self.progress_bars:
            progress.setValue(0)
            status.setText("Ready")
            status.setStyleSheet("")
        
        self.console.clear()
        self.clear_results_tab()
        
        # Disable UI during processing
        self.setEnabled(False)
        self.start_btn.setText("Processing...")
        QApplication.processEvents()
        
        # Execute first step
        self.execute_step(0, self.generate_smiles, formula, charge, "nSMILES.txt")
    
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
        
        # Mark step as completed
        self.current_step_results[step] = True
        self.update_progress(step, 100, "Completed successfully!")
        QApplication.processEvents()
        
        # Get current parameters
        formula = self.formula_input.text().strip()
        charge = self.charge_input.value()
        precursor_features = [self.precursor_list.item(i).text() 
                            for i in range(self.precursor_list.count())]
        branched = self.branched_check.isChecked()
        ring = self.ring_check.isChecked()
        target_mass = self.mass_input.value()
        tolerance = self.tolerance_input.value()
        
        # Execute next step
        if step == 0:  # Neutral SMILES generated
            self.execute_step(1, self.filter_smiles, "nSMILES.txt", formula, charge, 
                            precursor_features, branched, ring, "ParentRelatedSMILES.txt")
        elif step == 1:  # Filtered by precursor
            self.execute_step(2, self.generate_charged_smiles, formula, charge, 
                            "ParentRelatedSMILES.txt", "chargedSMILES.txt")
        elif step == 2:  # Charged SMILES generated
            self.execute_step(3, self.filter_charged_smiles_by_mass, formula, charge, 
                            "chargedSMILES.txt", target_mass, tolerance, "filteredchargedSMILES.txt")
        elif step == 3:  # Filtered by mass
            self.execute_step(4, self.generate_and_display_molecules, formula, charge, 
                            "filteredchargedSMILES.txt")
        elif step == 4:  # All done
            self.log_message("\nProcessing completed successfully!")
            self.cleanup_after_processing(success=True)
    
    def generate_and_display_molecules(self, formula, charge, input_file, progress_callback=None):
        try:
            input_path = os.path.join(self.output_dir, input_file)
            self.molecules, self.data = self.smiles_to_molecules(formula, charge, input_path, progress_callback)
            
            # Force GUI update
            QApplication.processEvents()
            self.display_results()
            
            return True
        except Exception as e:
            if progress_callback:
                progress_callback(100, f"Error: {str(e)}")
            return False
    
    def display_results(self):
        self.clear_results_tab()
        
        if not self.data:
            no_results = QLabel("No molecules were generated that match all criteria")
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
            
            # Display molecule info
            info_layout = QVBoxLayout()
            info_layout.addWidget(QLabel(f"<b>Formula:</b> {formula}"))
            info_layout.addWidget(QLabel(f"<b>Mass:</b> {mass}"))
            info_layout.addWidget(QLabel(f"<b>SMILES:</b> {smiles}"))
            info_layout.addStretch()
            layout.addLayout(info_layout)
            
            group.setLayout(layout)
            self.results_layout.addWidget(group)
        
        # Ensure scroll area updates
        self.results_area.verticalScrollBar().setValue(0)
        QApplication.processEvents()
    
    def clear_results_tab(self):
        # Properly clear previous results
        while self.results_layout.count():
            item = self.results_layout.takeAt(0)
            widget = item.widget()
            if widget:
                widget.deleteLater()
        QApplication.processEvents()
    
    def cleanup_after_processing(self, success=False):
        # Re-enable UI
        self.setEnabled(True)
        self.start_btn.setText("Start Processing")
        
        if success:
            # Ensure all steps show as completed
            for step in range(5):
                if step not in self.current_step_results:
                    self.update_progress(step, 100, "Completed")
            
            # Switch to results tab and force update
            self.central_widget.setCurrentIndex(2)
            QApplication.processEvents()
    
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

    # Chemistry functions
    def generate_smiles(self, FORMULA, CHARGE, output_file=None, progress_callback=None):
        output_dir = f"OutputFiles_{FORMULA}_Charge_{CHARGE}"
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, output_file) if output_file else None

        atoms = self.parse_formula(FORMULA)
        graphs = self.generate_graphs_lazy(atoms, progress_callback)
        smiles_set = set()
        total = 0

        for g in graphs:
            mol = self.graph_to_rdkit_mol(g)
            try:
                smiles = Chem.MolToSmiles(mol, canonical=True)
                smiles_set.add(smiles)
                total += 1
                if progress_callback and total % 10 == 0:
                    progress = 50 + int((total % 50))
                    progress_callback(progress, f"Generated {len(smiles_set)} SMILES")
            except:
                continue

        if output_path:
            with open(output_path, "w") as f:
                for smiles in sorted(smiles_set):
                    f.write(smiles + "\n")

        return len(smiles_set)

    def parse_formula(self, FORMULA):
        elements = re.findall(r'([A-Z][a-z]*)(\d*)', FORMULA)
        atoms = []
        for elem, count in elements:
            count = int(count) if count else 1
            atoms.extend([elem] * count)
        return atoms

    def generate_graphs_lazy(self, atoms, progress_callback=None):
        G_base = self.atoms_graph(atoms)
        n = len(atoms)
        possible_edges = list(combinations(range(n), 2))
        total_combinations = math.comb(len(possible_edges), n - 1)
        processed = 0

        for edge_comb in combinations(possible_edges, n - 1):
            g = G_base.copy()
            valid = True

            for (i, j) in edge_comb:
                a1, a2 = g.nodes[i]['element'], g.nodes[j]['element']
                if (a1, a2) in self.bond_orders or (a2, a1) in self.bond_orders:
                    g.add_edge(i, j, order=1)
                else:
                    valid = False
                    break

            if valid and nx.is_tree(g) and self.is_valid(g):
                expanded = self.expand_bond_orders(g)
                for eg in expanded:
                    if self.is_valid(eg):
                        yield eg
            
            processed += 1
            if progress_callback and processed % 100 == 0:
                progress = min(100, int(processed / total_combinations * 50))
                progress_callback(progress, f"Generating graphs: {processed}/{total_combinations}")

    def atoms_graph(self, atoms):
        G = nx.Graph()
        for i, atom in enumerate(atoms):
            G.add_node(i, element=atom)
        return G

    def is_valid(self, graph):
        for node in graph.nodes:
            atom = graph.nodes[node]['element']
            valence = sum(data['order'] for _, _, data in graph.edges(node, data=True))
            if valence > self.max_valence[atom]:
                return False
        return True

    def expand_bond_orders(self, graph):
        graphs = []
        edges = list(graph.edges(data=True))
        bond_options = []
        for (i, j, data) in edges:
            a1 = graph.nodes[i]['element']
            a2 = graph.nodes[j]['element']
            allowed_orders = self.bond_orders.get((a1, a2)) or self.bond_orders.get((a2, a1))
            bond_options.append([(i, j, order) for order in allowed_orders])
        for bond_combination in product(*bond_options):
            g = graph.copy()
            for (i, j, order) in bond_combination:
                g[i][j]['order'] = order
            graphs.append(g)
        return graphs

    def graph_to_rdkit_mol(self, graph):
        rw_mol = Chem.RWMol()
        node_to_idx = {}
        for node in graph.nodes:
            atom = Chem.Atom(graph.nodes[node]['element'])
            idx = rw_mol.AddAtom(atom)
            node_to_idx[node] = idx
        for i, j, data in graph.edges(data=True):
            rw_mol.AddBond(node_to_idx[i], node_to_idx[j], self.number_to_bondtype(data['order']))
        return rw_mol

    def filter_smiles(self, input_file, FORMULA, CHARGE, PRECURSOR_FEATURES, branched, ring, output_file=None, progress_callback=None):
        output_dir = f"OutputFiles_{FORMULA}_Charge_{CHARGE}"
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, output_file) if output_file else None
        input_path = os.path.join(output_dir, input_file)

        with open(input_path, "r") as f:
            smiles_list = [line.strip() for line in f if line.strip()]

        filtered = []
        total = len(smiles_list)

        for i, smi in enumerate(smiles_list):
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                continue

            if not ring and mol.GetRingInfo().NumRings() > 0:
                continue

            if not branched and not self.is_linear(mol):
                continue

            if self.contains_feature(mol, PRECURSOR_FEATURES):
                filtered.append(smi)

            if progress_callback and i % 10 == 0:
                progress = min(100, int(i / total * 100))
                progress_callback(progress, f"Filtered {i}/{total} SMILES")

        if output_path:
            with open(output_path, "w") as f:
                for smi in filtered:
                    f.write(smi + "\n")

        return len(filtered)

    def is_linear(self, mol):
        degrees = [atom.GetDegree() for atom in mol.GetAtoms()]
        num_ends = degrees.count(1)
        num_middle = degrees.count(2)
        return num_ends == 2 and (num_ends + num_middle == len(degrees))

    def contains_feature(self, mol, PRECURSOR_FEATURES):
        for feature in PRECURSOR_FEATURES:
            pattern = Chem.MolFromSmarts(feature)
            if pattern and mol.HasSubstructMatch(pattern):
                return True
        return False

    def generate_charged_smiles(self, FORMULA, CHARGE, input_file, output_file=None, progress_callback=None):
        output_dir = f"OutputFiles_{FORMULA}_Charge_{CHARGE}"
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, output_file) if output_file else None
        input_path = os.path.join(output_dir, input_file)

        with open(input_path, "r") as f:
            smiles_list = [line.strip() for line in f if line.strip()]

        charged_smiles = []
        total = len(smiles_list)

        for i, smi in enumerate(smiles_list):
            mol = Chem.MolFromSmiles(smi, sanitize=False)
            if mol is None:
                continue
            for atom in mol.GetAtoms():
                if atom.GetSymbol() in ["C", "O"]:
                    mol_copy = Chem.RWMol(mol)
                    atom_copy = mol_copy.GetAtomWithIdx(atom.GetIdx())
                    atom_copy.SetFormalCharge(CHARGE)
                    try:
                        charged_smi = Chem.MolToSmiles(mol_copy, canonical=True)
                        charged_smiles.append(charged_smi)
                    except:
                        continue

            if progress_callback and i % 10 == 0:
                progress = min(100, int(i / total * 100))
                progress_callback(progress, f"Generated {len(charged_smiles)} charged SMILES")

        if output_path:
            with open(output_path, "w") as f:
                for smi in charged_smiles:
                    f.write(smi + "\n")

        return len(charged_smiles)

    def filter_charged_smiles_by_mass(self, FORMULA, CHARGE, input_file, TARGET_MASS, tolerance, output_file=None, progress_callback=None):
        output_dir = f"OutputFiles_{FORMULA}_Charge_{CHARGE}"
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, output_file) if output_file else None
        input_path = os.path.join(output_dir, input_file)

        with open(input_path, "r") as f:
            smiles_list = [line.strip() for line in f if line.strip()]

        filtered_charged = []
        total = len(smiles_list)

        for i, smi in enumerate(smiles_list):
            mol = Chem.MolFromSmiles(smi)
            if mol:
                mass = Descriptors.ExactMolWt(mol)
                if abs(mass - TARGET_MASS) <= tolerance:
                    filtered_charged.append(smi)

            if progress_callback and i % 10 == 0:
                progress = min(100, int(i / total * 100))
                progress_callback(progress, f"Filtered {i}/{total} by mass")

        if output_path:
            with open(output_path, "w") as f:
                for smi in filtered_charged:
                    f.write(smi + "\n")

        return len(filtered_charged)

    def smiles_to_molecules(self, FORMULA, CHARGE, input_file, progress_callback=None):
        output_dir = f"OutputFiles_{FORMULA}_Charge_{CHARGE}"
        os.makedirs(output_dir, exist_ok=True)
        input_path = os.path.join(output_dir, input_file)

        with open(input_path, "r") as f:
            smiles_list = [line.strip() for line in f if line.strip()]

        molecules = []
        data = []
        seen_canonical_smiles = set()
        total = len(smiles_list)

        for i, smiles in enumerate(smiles_list):
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                canonical = Chem.MolToSmiles(mol, canonical=True)
                if canonical in seen_canonical_smiles:
                    continue
                seen_canonical_smiles.add(canonical)

                mass = f"{Descriptors.ExactMolWt(mol):.4f}"
                formula = rdMolDescriptors.CalcMolFormula(mol)
                data.append((smiles, mass, formula))
                molecules.append(mol)

                # Generate molecule image
                mol_img = Draw.MolToImage(mol, size=(300, 200))

                # Create annotated image
                annotated_img = self.create_annotated_image(mol_img, formula, mass, smiles)
                
                # Save image
                safe_smiles = self.sanitize_filename(smiles)
                filename = f"mol_{i + 1}_{safe_smiles}.png"
                output_path = os.path.join(output_dir, filename)
                annotated_img.save(output_path)

            if progress_callback and i % 5 == 0:
                progress = min(100, int(i / total * 100))
                progress_callback(progress, f"Processed {i}/{total} molecules")

        return molecules, data

    def create_annotated_image(self, mol_img, formula, mass, smiles):
        IMG_SIZE = (300, 200)
        ANNOTATION_HEIGHT = 60
        FONT_SIZE = 14

        total_height = IMG_SIZE[1] + ANNOTATION_HEIGHT
        annotated_img = Image.new("RGB", (IMG_SIZE[0], total_height), "white")
        annotated_img.paste(mol_img, (0, 0))

        draw = ImageDraw.Draw(annotated_img)
        try:
            font = ImageFont.truetype("arial.ttf", FONT_SIZE)
        except:
            font = ImageFont.load_default()

        annotation1 = f"{formula} | Mass: {mass}"
        annotation2 = smiles

        text_y = IMG_SIZE[1] + 5
        for text in [annotation1, annotation2]:
            bbox = draw.textbbox((0, 0), text, font=font)
            text_width = bbox[2] - bbox[0]
            x_text = (IMG_SIZE[0] - text_width) // 2
            draw.text((x_text, text_y), text, fill="black", font=font)
            text_y += FONT_SIZE + 2

        return annotated_img

    def number_to_bondtype(self, order):
        if order == 1: return Chem.BondType.SINGLE
        if order == 2: return Chem.BondType.DOUBLE
        if order == 3: return Chem.BondType.TRIPLE
        raise ValueError("Invalid bond order")

if __name__ == "__main__":
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    window = ForMileSGUI()
    window.show()
    sys.exit(app.exec_())