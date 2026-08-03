# Python-cellphonedb-workflow-template

CellPhoneDB、物种转换及通讯结果可视化代码素材模板。当前版本以保留原分析实现为主，尚未合并成一个统一正式入口。

## 主要代码

- `cellphonedb_run_mouse.py`：小鼠数据 CellPhoneDB 命令行流程。
- `mouse2human.py`：小鼠到人同源基因转换。
- `5_bubble_group.py`、`Chord_diagram.py`：结果可视化。
- `versions/`：来源不同的 CellPhoneDB 运行版本。
- `examples/`：原项目命令示例。
- `3_cellphonedb.py`：原 Streamlit 平台代码，仅用于参考。

## 使用前注意

历史脚本仍包含固定数据库路径、平台依赖及删除既有输出目录的逻辑。新项目使用时应先复制需要的代码到项目 `workflow/`，配置数据库和输入路径，并确认输出目录后再运行。

原始副本保存在 `../R-Python-bioinformatics-code-archive/legacy/original/05_cellphonedb/`。

## 其他人员来源代码

`source_catalog/` 收录 28 份 CellPhoneDB 语义副本，`SOURCE_CATALOG.tsv` 记录分析场景、owner、原文件名、来源路径和 SHA-256。根目录精选脚本用于快速选代码；来源目录用于比较不同项目实现，不建议不检查路径和参数就直接运行。
