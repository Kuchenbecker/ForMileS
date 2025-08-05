import subprocess
import webbrowser
import time

# Abre o navegador na interface
url = "http://localhost:8501"
webbrowser.open(url)

# Executa o app Streamlit
subprocess.run(["streamlit", "run", "app.py"], shell=True)
