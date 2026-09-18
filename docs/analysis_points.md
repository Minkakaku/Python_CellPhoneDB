# CellPhoneDB 分析点说明

## 01 物种与输入准备

CellPhoneDB 需要细胞条形码、细胞类型标签和与数据库标识一致的表达矩阵。人数据通常使用 HGNC symbol；小鼠数据先通过人工审核的同源表转换为人基因。本模板对多对一同源基因求和，避免旧代码只保留第一条映射造成表达丢失。

## 02 统计通讯

正式入口使用 CellPhoneDB v5 `cpdb_statistical_analysis_method.call()`：

1. 计算多亚基配体/受体复合物的表达；
2. 按细胞类型标签置换建立 null distribution；
3. 输出 means、pvalues、significant_means 和 deconvoluted 表；
4. 可用 microenvironment 文件限制空间上可能相遇的细胞群。

置换单位是细胞标签，因此结果用于候选通讯筛选，不等同于跨样本因果或组间推断。

## 03 可视化

`ktplotspy` 根据稳定的 means/pvalues 表绘制：

- 细胞类型之间显著互作数量热图；
- 配体-受体气泡图。

历史 chord diagram 和 Streamlit 页面保留在 `legacy/original/`。Chord 图是展示方式，不应改变统计筛选标准。
