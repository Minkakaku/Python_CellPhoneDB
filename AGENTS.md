# CellPhoneDB 项目约定

- 本仓库是独立分析点，正式入口仅为根目录 cellphonedb_workflow.py。
- 输入为主流程产生的注释后 h5ad；项目路径和参数通过命令行传入。
- 小鼠输入必须提供同源映射，many-to-one 表达求和。
- 最终图成对保存 PNG/PDF，结果表和 README_results_CN.md 放同一输出目录。
- CellPhoneDB 标签置换 P 值不能解释为跨生物学样本显著性。
- legacy/original/ 和 source_catalog/ 保持原字节。