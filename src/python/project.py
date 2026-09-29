"""项目配置驱动已有公共 CLI；dry-run 不加载分析依赖。"""
from pathlib import Path
import contextlib
import json
import runpy
import sys
import yaml

def merge_config(base, override):
    result = dict(base)
    for key, value in override.items():
        result[key] = merge_config(result[key], value) if isinstance(value, dict) and isinstance(result.get(key), dict) else value
    return result

def load_method(project_dir, method):
    project = Path(project_dir).resolve()
    workspace = project.parents[1]
    config = {}
    for path in (workspace / "config/default.yaml", workspace / "config/reference.yaml", project / "config/config.yaml"):
        config = merge_config(config, yaml.safe_load(path.read_text()) or {})
    options = config["methods"][method]
    options = dict(options)
    path_keys = {"input", "input-h5ad", "output", "output-dir", "object-out", "object-dir", "gtf", "gmt", "tf-list", "databases", "motif-annotations", "database", "homology", "microenvs"}
    for key in path_keys & options.keys():
        value = options[key]
        if value is not None:
            resolve = lambda item: str((project / Path(item).expanduser()).resolve())
            options[key] = [resolve(item) for item in value] if isinstance(value, list) else resolve(value)
    args = []
    for key, value in options.items():
        if value is None or value is False:
            continue
        args.append("--" + key)
        if value is not True:
            args.extend(str(v) for v in (value if isinstance(value, list) else [value]))
    return workspace, project, args

def run_project_method(project_dir, method):
    workspace, project, args = load_method(project_dir, method)
    args.extend(x for x in sys.argv[1:] if x != "--dry-run")
    script = workspace / "src/python" / (method + "_workflow.py")
    if "--dry-run" in sys.argv:
        print(json.dumps({"cwd": str(project), "script": str(script), "args": args}, ensure_ascii=False))
        return
    import os
    from datetime import datetime
    (project / "logs").mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S_%f")
    previous_cwd, previous_argv = Path.cwd(), sys.argv
    try:
        os.chdir(project)
        sys.argv = [str(script), *args]
        with (project / "logs" / f"{method}_{stamp}.log").open("w") as log:
            log.write(json.dumps({"cwd": str(project), "args": args}, ensure_ascii=False) + "\n")
            log.flush()
            with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log):
                runpy.run_path(str(script), run_name="__main__")
    finally:
        os.chdir(previous_cwd)
        sys.argv = previous_argv
