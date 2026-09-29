from pathlib import Path
import sys

project = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(project.parents[1] / "src/python"))
from project import run_project_method

run_project_method(project, "cellphonedb")
