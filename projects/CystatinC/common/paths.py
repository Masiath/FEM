"""output location for every script in this project."""
from pathlib import Path
OUT = Path(__file__).resolve().parent / "results"
OUT.mkdir(exist_ok=True)
OUT = str(OUT)
