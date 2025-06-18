from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import Descriptors
from PIL import Image, ImageDraw, ImageFont
import os
import re

IMG_SIZE = (300, 200)
ANNOTATION_HEIGHT = 60
FONT_SIZE = 14

def sanitize_filename(smiles):
    return re.sub(r'[^a-zA-Z0-9._-]', '_', smiles)

def smiles_to_molecules(FORMULA, CHARGE, input_file, as_svg=False):
    output_dir = f"OutputFiles_{FORMULA}_Charge_{CHARGE}"
    input_path = os.path.join(output_dir, input_file)
    os.makedirs(output_dir, exist_ok=True)

    with open(input_path, "r") as f:
        smiles_list = [line.strip() for line in f if line.strip()]

    molecules = []
    data = []
    image_paths = []
    seen_canonical_smiles = set()
    BATCH_SIZE = 100

    for idx in range(0, len(smiles_list), BATCH_SIZE):
        batch = smiles_list[idx:idx+BATCH_SIZE]
        for smi in batch:
            mol = Chem.MolFromSmiles(smi)
            if mol:
                canonical = Chem.MolToSmiles(mol, canonical=True)
                if canonical in seen_canonical_smiles:
                    continue
                seen_canonical_smiles.add(canonical)

                mass = f"{Descriptors.ExactMolWt(mol):.4f}"
                formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
                data.append((smi, mass, formula))
                molecules.append(mol)

                safe_smiles = sanitize_filename(smi)
                png_path = os.path.join(output_dir, f"mol_{idx + 1}_{safe_smiles}.png")
                svg_path = os.path.join(output_dir, f"mol_{idx + 1}_{safe_smiles}.svg")

                # Generate PNG
                mol_img = Draw.MolToImage(mol, size=IMG_SIZE)
                total_height = IMG_SIZE[1] + ANNOTATION_HEIGHT
                annotated_img = Image.new("RGB", (IMG_SIZE[0], total_height), "white")
                annotated_img.paste(mol_img, (0, 0))

                draw = ImageDraw.Draw(annotated_img)
                try:
                    font = ImageFont.truetype("arial.ttf", FONT_SIZE)
                except:
                    font = ImageFont.load_default()

                annotation1 = f"{formula} | Mass: {mass}"
                annotation2 = smi
                text_y = IMG_SIZE[1] + 5
                for text in [annotation1, annotation2]:
                    bbox = draw.textbbox((0, 0), text, font=font)
                    text_width = bbox[2] - bbox[0]
                    x_text = (IMG_SIZE[0] - text_width) // 2
                    draw.text((x_text, text_y), text, fill="black", font=font)
                    text_y += FONT_SIZE + 2

                annotated_img.save(png_path)
                image_paths.append(png_path)

                # Generate SVG if requested
                if as_svg:
                    drawer = Draw.MolDraw2DSVG(IMG_SIZE[0], IMG_SIZE[1])
                    drawer.DrawMolecule(mol)
                    drawer.FinishDrawing()
                    svg = drawer.GetDrawingText()
                    with open(svg_path, "w") as f:
                        f.write(svg)

    if not molecules:
        print("No valid molecules to display.")

    return molecules, data, image_paths