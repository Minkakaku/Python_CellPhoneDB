#!/usr/bin/env python
"""CellPhoneDB v5 独立分析入口。

输入是主单细胞流程产生的 ``.h5ad``，输出 CellPhoneDB 统计通讯结果、
热图、气泡图和中文结果说明。小鼠数据会先按同源表转换成人基因；
多个小鼠基因映射到同一人基因时表达量相加。
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
from scipy import sparse


NPG_SEQUENTIAL = ["#FFF5F0", "#F39B7F", "#E64B35", "#7A0019"]


# 功能：读取主流程输出，并选取 CellPhoneDB 使用的表达层。
def read_input(path: Path, celltype_col: str, layer: str | None, use_raw: bool) -> ad.AnnData:
    adata = ad.read_h5ad(path)
    if use_raw:
        if adata.raw is None:
            raise ValueError("--use-raw 已启用，但 h5ad 中没有 adata.raw。")
        raw = adata.raw.to_adata()
        raw.obs = adata.obs.copy()
        adata = raw
    elif layer:
        adata.X = adata.layers[layer].copy()
    if celltype_col not in adata.obs:
        raise KeyError(f"obs 中缺少细胞类型列：{celltype_col}")
    adata.obs_names_make_unique()
    adata.var_names_make_unique()
    adata.obs[celltype_col] = adata.obs[celltype_col].astype("category")
    return adata


# 功能：把多个小鼠基因映射到同一人基因时按表达量求和。
def convert_mouse_to_human(
    adata: ad.AnnData,
    mapping_path: Path,
    mouse_col: str,
    human_col: str,
) -> ad.AnnData:
    mapping = pd.read_csv(mapping_path, sep=None, engine="python")
    mapping = mapping[[mouse_col, human_col]].dropna().drop_duplicates()
    mapping = mapping[mapping[mouse_col].isin(adata.var_names)].copy()
    if mapping.empty:
        raise ValueError("同源映射表与 h5ad 基因没有交集。")

    source_index = pd.Series(np.arange(adata.n_vars), index=adata.var_names)
    human_genes = pd.Index(pd.unique(mapping[human_col].astype(str)))
    human_index = pd.Series(np.arange(len(human_genes)), index=human_genes)
    rows = source_index.loc[mapping[mouse_col]].to_numpy()
    cols = human_index.loc[mapping[human_col].astype(str)].to_numpy()

    # 映射矩阵完成 many-to-one 求和，避免只保留第一个同源基因。
    gene_map = sparse.csr_matrix(
        (np.ones(len(rows)), (rows, cols)),
        shape=(adata.n_vars, len(human_genes)),
    )
    expression = adata.X @ gene_map
    if not sparse.issparse(adata.X):
        expression = np.asarray(expression)

    converted = ad.AnnData(
        X=expression,
        obs=adata.obs.copy(),
        var=pd.DataFrame(index=human_genes),
    )
    converted.uns["homology_conversion"] = {
        "source_species": "mouse",
        "target_species": "human",
        "mapping_file": str(mapping_path),
        "mapped_mouse_genes": int(mapping[mouse_col].nunique()),
        "human_genes": int(len(human_genes)),
        "aggregation": "sum_many_to_one",
    }
    return converted


# 功能：写出 CellPhoneDB 要求的 cell/cell_type 两列 metadata。
def write_metadata(adata: ad.AnnData, celltype_col: str, path: Path) -> None:
    pd.DataFrame(
        {"Cell": adata.obs_names, "cell_type": adata.obs[celltype_col].astype(str)}
    ).to_csv(path, sep="\t", index=False)


# 功能：运行官方 statistical_analysis；返回值按 v5 API 作为字典处理。
def run_cellphonedb(args: argparse.Namespace, counts_path: Path, metadata_path: Path) -> dict:
    from cellphonedb.src.core.methods import cpdb_statistical_analysis_method

    # 置换检验回答“某细胞对是否比其余细胞对更特异”，不是跨样本条件检验。
    return cpdb_statistical_analysis_method.call(
        cpdb_file_path=str(args.database),
        meta_file_path=str(metadata_path),
        counts_file_path=str(counts_path),
        counts_data=args.counts_data,
        microenvs_file_path=None if args.microenvs is None else str(args.microenvs),
        iterations=args.iterations,
        threshold=args.threshold,
        threads=args.threads,
        debug_seed=args.seed,
        result_precision=4,
        pvalue=args.pvalue,
        subsampling=args.subsampling,
        subsampling_log=False,
        subsampling_num_cells=args.subsampling_num_cells,
        separator="|",
        debug=False,
        output_path=str(args.output),
        output_suffix="result",
    )


# 功能：把当前 Matplotlib 图成对保存为 PDF 和 600 dpi PNG。
def save_matplotlib(stem: Path, width: float, height: float) -> None:
    figure = plt.gcf()
    figure.set_size_inches(width, height)
    figure.tight_layout()
    figure.savefig(stem.with_suffix(".pdf"), bbox_inches="tight")
    figure.savefig(stem.with_suffix(".png"), dpi=600, bbox_inches="tight")
    plt.close(figure)


# 功能：使用 ktplotspy 绘制通讯数量热图和全局配体-受体气泡图。
def plot_results(
    adata: ad.AnnData,
    celltype_col: str,
    means: pd.DataFrame,
    pvalues: pd.DataFrame,
    output: Path,
) -> None:
    import ktplotspy as kpy

    kpy.plot_cpdb_heatmap(
        adata=adata,
        pvals=pvalues,
        celltype_key=celltype_col,
        figsize=(10, 9),
        title="",
        symmetrical=False,
    )
    cmap = LinearSegmentedColormap.from_list("npg_sequential", NPG_SEQUENTIAL)
    for axis in plt.gcf().axes:
        for artist in [*axis.collections, *axis.images]:
            artist.set_cmap(cmap)
    save_matplotlib(output / "interaction_count_heatmap", 10, 9)

    bubble = kpy.plot_cpdb(
        adata=adata,
        cell_type1=".",
        cell_type2=".",
        means=means,
        pvals=pvalues,
        celltype_key=celltype_col,
        figsize=(14, 10),
        title="",
        max_size=6,
        highlight_size=0.75,
        standard_scale=True,
    )
    from plotnine import element_blank, scale_color_gradientn, theme

    # 通讯气泡点图不显示主/次内部网格，保留原有外框。
    bubble += scale_color_gradientn(colors=NPG_SEQUENTIAL) + theme(
        panel_grid_major=element_blank(), panel_grid_minor=element_blank()
    )
    bubble.save(str(output / "ligand_receptor_bubble.pdf"), width=14, height=10, verbose=False)
    bubble.save(
        str(output / "ligand_receptor_bubble.png"),
        width=14,
        height=10,
        dpi=600,
        verbose=False,
    )


# 功能：记录输入、参数、统计含义和输出限制，便于项目复核。
def write_result_notes(args: argparse.Namespace, counts_path: Path, tables: list[str]) -> None:
    (args.output / "README_results_CN.md").write_text(
        "\n".join(
            [
                "# CellPhoneDB 结果说明",
                "",
                f"- 输入 h5ad：{args.input}",
                f"- 实际分析对象：{counts_path}",
                f"- 数据库：{args.database}",
                f"- 细胞类型列：{args.celltype_col}",
                f"- species：{args.species}；gene identifier：{args.counts_data}",
                f"- iterations：{args.iterations}；threshold：{args.threshold}；nominal P：{args.pvalue}",
                f"- 输出表：{', '.join(tables)}",
                "- statistical_analysis 的 P 值来自细胞标签置换，用于判断细胞对特异性。",
                "- 它不是生物学重复层面的 CON/EXP 差异检验；条件比较应按样本分别运行并在样本层面汇总。",
                "- 多亚基复合物要求各亚基达到表达阈值，互作均值由最低表达亚基限制。",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


# 功能：定义独立脚本的输入对象、数据库、物种和置换参数。
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run CellPhoneDB v5 from an annotated h5ad.")
    parser.add_argument("--input", type=Path, required=True, help="主流程输出的 h5ad")
    parser.add_argument("--database", type=Path, required=True, help="cellphonedb.zip")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--celltype-col", default="cell_type")
    parser.add_argument("--species", choices=["human", "mouse"], default="human")
    parser.add_argument("--counts-data", default="hgnc_symbol")
    parser.add_argument("--layer", default=None, help="可选表达层；不要使用中心化 scale 数据")
    parser.add_argument("--use-raw", action="store_true", help="使用 adata.raw")
    parser.add_argument("--homology", type=Path, default=None)
    parser.add_argument("--mouse-col", default="mouse_gene")
    parser.add_argument("--human-col", default="human_gene")
    parser.add_argument("--microenvs", type=Path, default=None)
    parser.add_argument("--iterations", type=int, default=1000)
    parser.add_argument("--threshold", type=float, default=0.1)
    parser.add_argument("--pvalue", type=float, default=0.05)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--seed", type=int, default=1001)
    parser.add_argument("--subsampling", action="store_true")
    parser.add_argument("--subsampling-num-cells", type=int, default=15000)
    parser.add_argument("--skip-plots", action="store_true")
    return parser


# 功能：串联输入准备、同源转换、CellPhoneDB、绘图和结果说明。
def main() -> None:
    args = build_parser().parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    adata = read_input(args.input, args.celltype_col, args.layer, args.use_raw)

    if args.species == "mouse":
        if args.homology is None:
            raise ValueError("小鼠数据必须通过 --homology 提供 mouse-human 同源表。")
        adata = convert_mouse_to_human(adata, args.homology, args.mouse_col, args.human_col)

    # CellPhoneDB 直接读取 h5ad；保存实际传入的矩阵以固定输入层和同源转换结果。
    counts_path = args.output / "cellphonedb_input.h5ad"
    adata.write_h5ad(counts_path)
    metadata_path = args.output / "cell_metadata.tsv"
    write_metadata(adata, args.celltype_col, metadata_path)

    results = run_cellphonedb(args, counts_path, metadata_path)
    tables = []
    for name, value in results.items():
        if isinstance(value, pd.DataFrame):
            value.to_csv(args.output / f"{name}.csv", index=False)
            tables.append(name)

    if not args.skip_plots:
        plot_results(adata, args.celltype_col, results["means"], results["pvalues"], args.output)
    write_result_notes(args, counts_path, tables)


if __name__ == "__main__":
    main()
