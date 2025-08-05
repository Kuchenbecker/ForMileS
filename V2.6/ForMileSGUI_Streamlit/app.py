import streamlit as st
import json
import subprocess

CONFIG_PATH = "config.json"

st.set_page_config(page_title="ForMileS Generator", layout="centered")
st.title("🧪 ForMileS: Molecular Structure Generator")

# Campos do usuário
formula = st.text_input("Fórmula molecular (ex: C6H6)", "C7Cl")
scaffold = st.text_input("SMILES do precursor (ex: C1=CC=CC=C1)", "C1C=CC=CC=1")
charge = st.number_input("Carga (ex: -1, 0, +1)", value=-1, step=1)
mass = st.number_input("Massa alvo (Da)", value=125.016)
tol = st.number_input("Tolerância de massa (Da)", value=0.5)

st.markdown("### ⚙️ Opções estruturais")
branching = st.checkbox("Permitir ramificações", value=True)
cycles = st.checkbox("Permitir ciclos", value=True)
double_bonds = st.checkbox("Permitir duplas ligações", value=True)
triple_bonds = st.checkbox("Permitir triplas ligações", value=False)

max_dbl = st.number_input("Máximo de duplas ligações", value=3)
max_trpl = st.number_input("Máximo de triplas ligações", value=2)

st.markdown("### 💾 Saída")
save_xyz = st.checkbox("Salvar arquivos .xyz", value=True)
save_mol = st.checkbox("Salvar arquivos .mol", value=True)
save_svg = st.checkbox("Salvar imagens .svg", value=True)

if st.button("🚀 Executar ForMileS"):
    config = {
        "FORMULA": formula,
        "MOL_SCAFFOLD": [scaffold],
        "CHARGE": charge,
        "TARGET_MASS": mass,
        "TOLERANCE": tol,
        "PARAM_FILE": "parameters.json",
        "SAVE_XYZ": save_xyz,
        "SAVE_MOL": save_mol,
        "SAVE_SVG": save_svg,
        "ALLOW_BRANCHING": branching,
        "ALLOW_CYCLES": cycles,
        "ALLOW_DOUBLE_BONDS": double_bonds,
        "ALLOW_TRIPLE_BONDS": triple_bonds,
        "MAX_DOUBLE_BONDS": max_dbl,
        "MAX_TRIPLE_BONDS": max_trpl
    }

    # Salva config
    with open(CONFIG_PATH, "w") as f:
        json.dump(config, f, indent=2)

    with st.spinner("Executando geração molecular..."):
        result = subprocess.run(["python", "main.py"], capture_output=True, text=True)

    st.success("Execução concluída!")
    st.text_area("Log do programa:", result.stdout + "\n" + result.stderr, height=300)

    st.markdown("---")
    st.markdown("Arquivos gerados estão na pasta de saída correspondente ao nome da fórmula.")

