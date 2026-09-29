"""Run with python tests/test_project_layout.py; requires only PyYAML."""
from pathlib import Path
import runpy
import tempfile
import shutil
import yaml

ROOT = Path(__file__).resolve().parents[1]
api = runpy.run_path(str(ROOT / "src/python/project.py"))
with tempfile.TemporaryDirectory() as tmp:
    workspace = Path(tmp).resolve() / "workspace"
    shutil.copytree(ROOT / "config", workspace / "config")
    for name in ("a", "b"):
        shutil.copytree(ROOT / "projects/project1/config", workspace / "projects" / name / "config")
    a, b = workspace / "projects/a", workspace / "projects/b"
    cfg = yaml.safe_load((a / "config/config.yaml").read_text())
    method = next(iter(cfg["methods"]))
    cfg["methods"][method]["custom-test"] = "with spaces"
    (a / "config/config.yaml").write_text(yaml.safe_dump(cfg))
    _, first, args = api["load_method"](a, method)
    _, second, other = api["load_method"](b, method)
    assert first == a and second == b
    assert args[-2:] == ["--custom-test", "with spaces"] or args[args.index("--custom-test") + 1] == "with spaces"
    assert "--custom-test" not in other
    assert all(isinstance(value, str) for value in args)
    input_flag = "--input-h5ad" if "--input-h5ad" in args else "--input"
    assert Path(args[args.index(input_flag) + 1]).is_relative_to(a)
    # Execute a harmless stand-in to check cwd, argument passing, logs and restoration.
    script = workspace / "src/python" / (method + "_workflow.py")
    script.parent.mkdir(parents=True)
    script.write_text("from pathlib import Path\nPath('probe.txt').write_text(str(Path.cwd()))\n")
    previous = Path.cwd()
    api["run_project_method"](a, method)
    assert Path.cwd() == previous
    assert (a / "probe.txt").read_text() == str(a)
    assert list((a / "logs").glob("*.log"))
    assert api["merge_config"]({"a": {"x": 1, "y": [1, 2]}}, {"a": {"y": [3]}}) == {"a": {"x": 1, "y": [3]}}
print("Project isolation and configuration override checks passed")
