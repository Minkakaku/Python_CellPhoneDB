# Python-cellphonedb-workflow-template

独立 CellPhoneDB v5 分析点。正式代码只有根目录的 cellphonedb_workflow.py，直接接收主 Scanpy 流程生成的注释后 h5ad。

## 正式流程

1. 从 h5ad 选择非负、未中心化的表达层。
2. 小鼠数据按显式同源表转换成人基因；many-to-one 映射按表达量求和。
3. 生成 CellPhoneDB metadata 并运行 statistical_analysis。
4. 保存完整结果表、通讯数量热图、配体-受体气泡图和中文说明。

CellPhoneDB v5 的 cpdb_statistical_analysis_method.call() 返回结果字典；正式脚本按 means、pvalues、significant_means 等键保存，不能按四个位置返回值解包。

## 运行

    conda env create -f environment.yml
    conda activate cellphonedb-workflow
    pip install -e .

    cellphonedb-workflow --input annotated.h5ad --database cellphonedb.zip --output results/cellphonedb --celltype-col cell_type --species human

小鼠数据增加：

    cellphonedb-workflow --input annotated_mouse.h5ad --database cellphonedb.zip --output results/cellphonedb --celltype-col cell_type --species mouse --homology mouse_human_homology.tsv

映射表默认需要 mouse_gene 和 human_gene 两列。可用 --mouse-col、--human-col 修改列名。

## 关键解释

- 官方默认 threshold=0.1：配体/受体需在对应细胞群至少 10% 细胞表达。
- statistical analysis 的 P 值来自细胞标签置换，用于细胞对特异性，不是跨样本 CON/EXP 显著性。
- 条件比较应按生物学样本分别运行，再在样本层面汇总。
- 微环境文件可限制组织中不可能相遇的细胞类型组合。
- legacy/original/ 和 source_catalog/ 仅用于来源追溯，不是正式入口。

方法依据：CellPhoneDB 官方结果文档与方法论文（PMID: 32103204）。