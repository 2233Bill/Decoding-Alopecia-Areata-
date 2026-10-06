<a id="chinese"></a>

# 斑秃基因表达分析

**使用 GEO 基因表达数据，分析疾病相关表达差异，并比较患者与对照的分类模型。**

[English version ↓](#english)

## 项目概览

本项目研究基因表达特征能否区分患者与对照、哪些特征有助于分类，以及模型选择和正则化如何影响表现。完整分析在一份 [R Markdown 文件](final%20coding.Rmd)中串联了探索性分析、差异表达检验、特征筛选和分类建模。

**数据：** GEO 数据集 [GSE68801](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68801)，通过 `GEOquery` 读取 `GSE68801_series_matrix.txt.gz`。输入包括探针表达矩阵及疾病状态、年龄、性别元数据。代码将 `Normal` 编为 `Control`，其他疾病状态编为 `Patient`，因此目标是患者与对照的二分类。仓库未包含数据文件。

## 关键结果

| 分析 | 原文记录的结果 | 含义 |
| --- | --- | --- |
| [最终 LASSO 分类器](final%20coding.Rmd#L859-L987) | **AUC 0.981**；`alpha = 1`、`lambda = 0.001`；汇总 10 折交叉验证预测 | 最终比较中，LASSO 的表现最好。 |
| [初次 LASSO](final%20coding.Rmd#L184-L232) | 约 20% 留出集上的 **AUC 0.9748** | 选定表达特征在此次划分中具有较强区分能力。 |
| [500 个选定探针的 PCA](final%20coding.Rmd#L109-L132) | **PC1 解释 57.3% 的方差** | 描述利用疾病标签筛选后的特征中的变异。 |

**证据范围：** 上述数字来自原分析的文字记录，仓库未保存渲染结果或拟合模型，本次未独立复算指标。两个 AUC 的评估设置不同，不能据此认定模型获得了提升。部分验证步骤之外已完成特征筛选，分数可能偏乐观；PCA 也使用了依据标签选定的特征。详见[分析局限](#分析局限)。

## 分析流程

GEO 数据导入 → 标签与元数据整理 → 表达分布检查 → `limma` 差异表达检验 → LASSO 特征筛选 → 四类模型比较与调参 → 评估与特征解释。

## 关键分析与决策

1. **同时考虑统计显著性与效应大小。** 使用 `limma::lmFit()` 和 `eBayes()` 估计患者与对照的表达差异。候选探针同时满足 FDR 校正后的 `p < 0.05` 和 `|logFC| > 1`，结合证据强度与表达变化缩小特征范围，再通过火山图、MA 图和 Q–Q 图检查结果。

2. **先压缩特征，再比较模型。** 将表达矩阵转换为“样本 × 特征”，使用二项逻辑回归 LASSO，在训练部分通过 `cv.glmnet()` 选择 `lambda.min`，保留非零系数探针，再加入数值型年龄和编码后的性别供后续分类器使用。代码探索了这组联合输入，但没有单独验证人口学变量带来的增益。

3. **比较模型类型及参数敏感性。** 使用 `caret` 比较随机森林、正则化逻辑回归、径向核 SVM 和 kNN，以 10 折交叉验证和 ROC 指标调参。另通过 20 次 80/20 重复划分，以及 LASSO `lambda`、SVM cost 参数实验考察表现变化。最终 `glmnet` 模型固定 `alpha = 1`，采用纯 LASSO；这些实验的解释仍受验证局限约束。

4. **将预测与特征解释连接起来。** 通过 LASSO 系数、`caret::varImp()` 排名、特征列表交集和 `hgu133plus2.db` 探针注释探索候选特征。这些步骤提供进一步调查的线索，不能证明因果关系或确认生物标志物。

## 核心代码与图表

仓库目前只有分析源文件。以下图表由代码生成，尚无可直接嵌入的图片或 HTML 报告。

| 查看位置 | 展示内容 |
| --- | --- |
| [数据准备与探索性分析](final%20coding.Rmd#L29-L182) | 元数据编码、表达箱线图、差异表达、火山图、MA 图、Q–Q 图及 PCA |
| [LASSO 特征筛选](final%20coding.Rmd#L184-L242) | 分层划分、正则化选择、非零系数探针及联合输入 |
| [模型调参与最终比较](final%20coding.Rmd#L656-L987) | 参数网格、ROC 曲线与最终 AUC 比较图 |
| [评估与解释](final%20coding.Rmd#L990-L1197) | 系数、变量重要性图及指标计算 |
| [参数实验与注释](final%20coding.Rmd#L1249-L1437) | 重复留出实验的 AUC 分布与探针到基因符号的映射 |

## 工具与技术

| 用途 | 实际使用 |
| --- | --- |
| 数据读取与统计检验 | R、`GEOquery`、`limma` |
| 数据转换与可视化 | `dplyr`、`tidyr`、`reshape2`、`ggplot2`、基础 R |
| 特征选择与分类 | `glmnet`、`caret`、`randomForest`、`e1071`；`caret` 的 `svmRadial` 使用 `kernlab` 引擎 |
| ROC 评估与注释 | `pROC`、`hgu133plus2.db` |
| 分析记录 | R Markdown |

## 分析局限

- **验证中的信息泄漏：** 差异表达筛选在首次训练/测试划分前使用了全部样本，后续交叉验证和重复留出实验又复用固定的 LASSO 特征列表，未在各训练折内重新筛选。因此 AUC 仍属于探索性结果，不能证明对独立队列的泛化能力。
- **指标解释：** `Control` 是因子的第一个水平，`twoClassSummary` 的敏感度因此针对对照类别（[caret 文档](https://topepo.github.io/caret/measuring-performance.html)）。自定义 `QLIKE` 不等于二元对数损失。[混淆矩阵绘图](final%20coding.Rmd#L1199-L1227)使用手填矩阵覆盖计算结果，因此本 README 未将其数值作为已验证成果。
- **稳定性比较：** [留出实验与最终模型对比图](final%20coding.Rmd#L1159-L1185)重复填入同一个最终 AUC 20 次，并非 20 次独立评估，不能据此判断调参后模型的方差或稳定性。
- **特征排名：** 代码将 LASSO 系数绝对值与 `caret` 重要性分数取均值，但这些指标的尺度与含义不同，不能作为统一的生物学排名。SVM/kNN 的 `caret` 重要性使用独立于模型的筛选方法，与原文描述的置换过程不同（[caret 文档](https://topepo.github.io/caret/variable-importance.html)）。

## 仓库结构

```text
.
├── README.md          # 项目展示说明
└── final coding.Rmd   # 完整分析源文件
```

## 运行说明

**需要准备数据并核对原代码，本次未验证完整 HTML 渲染。**

1. 安装下列依赖。
2. 从上述 GEO 页面获取 `GSE68801_series_matrix.txt.gz`，将 `getGEO(filename = ...)` 中的 Windows 绝对路径改为本地路径。
3. 在 RStudio 中打开 `final coding.Rmd`，按顺序执行代码块。先核对性别元数据：初始代码使用 `gender:ch1 == "M"`，重复留出部分使用 `Sex:ch1 == "male"`，需依据实际数据核对字段及编码后再尝试完整运行。

<details>
<summary>依赖安装</summary>

```r
install.packages(c(
  "BiocManager", "R.utils", "reshape2", "ggplot2", "dplyr", "tidyr",
  "pheatmap", "RColorBrewer", "DT", "caret", "randomForest",
  "glmnet", "pROC", "e1071", "kernlab", "rmarkdown"
))
BiocManager::install(c("GEOquery", "limma", "hgu133plus2.db"))
```

以上涵盖原文件加载的包及 SVM 引擎依赖，仓库未锁定版本。另需注意：前面的调参文字记录写了 kNN `k = 9`，最终模型代码使用的是 `k = 3`。

</details>

## 项目体现的能力

原分析显示选定表达特征具有疾病相关信号，LASSO 在最终比较中取得最高的记录 AUC。项目展示了元数据整理、高维矩阵处理、统计筛选与正则化结合、多模型比较，以及通过诊断图表和特征解释沟通结果的能力；验证局限也说明了判断证据边界的重要性。

---

<a id="english"></a>

# Decoding Alopecia Areata

**An R analysis of disease-associated expression patterns and patient–control classification using GEO data.**

[中文版本 ↑](#chinese)

## Project Overview

This project asks whether gene-expression features can distinguish patients from controls, which features contribute to classification, and how model choice and regularization affect performance. It connects exploratory analysis, differential-expression testing, feature selection, and classification in one [R Markdown notebook](final%20coding.Rmd).

**Data:** GEO accession [GSE68801](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE68801), read from `GSE68801_series_matrix.txt.gz` with `GEOquery`. Inputs include a probe-expression matrix and disease-status, age, and gender metadata. The code maps `Normal` to `Control` and other disease-status values to `Patient`; it evaluates binary classification rather than disease subtypes. The dataset is not included in this repository.

## Key Results

| Analysis | Result recorded in the notebook | Interpretation |
| --- | --- | --- |
| [Final LASSO classifier](final%20coding.Rmd#L859-L987) | **AUC 0.981**; `alpha = 1`, `lambda = 0.001`; pooled predictions from 10-fold cross-validation | The final comparison reports LASSO as the strongest classifier. |
| [Initial LASSO](final%20coding.Rmd#L184-L232) | **AUC 0.9748** on the approximately 20% hold-out set | Selected expression features show strong discrimination in this split. |
| [PCA of 500 selected probes](final%20coding.Rmd#L109-L132) | **PC1 explains 57.3%** of variance | Describes variation within the features selected using disease labels. |

**Evidence boundary:** These numbers are recorded in the notebook's narrative; rendered outputs and fitted models are not committed, and the metrics have not been independently reproduced here. The two AUCs use different evaluation settings and do not establish an improvement of one model over the other. Feature selection occurs outside parts of the validation process, so the scores may be optimistic. The PCA is also based on label-selected features. See [limitations](#limitations).

## Analysis Approach

GEO import → label and metadata preparation → expression diagnostics → `limma` differential-expression testing → LASSO feature selection → four-model comparison and tuning → evaluation and feature interpretation.

## Key Analysis & Reasoning

1. **Combine statistical significance with effect size.** `limma::lmFit()` and `eBayes()` estimate patient–control expression differences. Candidate probes must satisfy FDR-adjusted `p < 0.05` **and** `|logFC| > 1`, narrowing the feature set using both evidence strength and expression change. Volcano, MA, and Q–Q plots inspect the resulting patterns.

2. **Reduce the feature set before broader model comparison.** The expression matrix is transposed into samples × predictors. Binomial LASSO uses `cv.glmnet()` to select `lambda.min` on the training partition and retains non-zero coefficients. These probes are then combined with numeric age and encoded gender for subsequent classifiers; no separate comparison establishes the benefit of the demographic variables.

3. **Compare model families and parameter sensitivity.** `caret` trains random forest, regularized logistic regression, radial SVM, and kNN with 10-fold cross-validation and ROC-based tuning. Separate experiments use 20 repeated 80/20 splits and vary LASSO `lambda` and SVM cost to inspect variation across splits and settings. The final `glmnet` model fixes `alpha = 1`, making it pure LASSO. These experiments remain subject to the validation limitations below.

4. **Connect prediction with feature interpretation.** The notebook examines LASSO coefficients, `caret::varImp()` rankings, overlaps between feature lists, and probe-to-symbol mapping with `hgu133plus2.db`. These provide candidates for investigation, without establishing causality or validated biomarkers.

## Code & Visualisations

The repository contains analysis source only. The plots below are generated by notebook chunks; no exported images or HTML report are available to embed.

| Where to look | What it shows |
| --- | --- |
| [Data preparation & EDA](final%20coding.Rmd#L29-L182) | Metadata recoding, expression boxplots, differential expression, volcano/MA/Q–Q plots, and PCA |
| [LASSO feature selection](final%20coding.Rmd#L184-L242) | Stratified split, regularization selection, non-zero probes, and combined predictors |
| [Model tuning & final comparison](final%20coding.Rmd#L656-L987) | Parameter grids, ROC curves, and final-model AUC chart |
| [Evaluation & interpretation](final%20coding.Rmd#L990-L1197) | Coefficient/importance plots and metric calculations |
| [Parameter experiments & annotation](final%20coding.Rmd#L1249-L1437) | Repeated hold-out AUC distributions and probe-to-gene mapping |

## Tools & Technologies

| Purpose | Used in the analysis |
| --- | --- |
| Data access & statistical testing | R, `GEOquery`, `limma` |
| Data transformation & plotting | `dplyr`, `tidyr`, `reshape2`, `ggplot2`, base R |
| Feature selection & classification | `glmnet`, `caret`, `randomForest`, `e1071`; `caret`'s `svmRadial` engine uses `kernlab` |
| ROC evaluation & annotation | `pROC`, `hgu133plus2.db` |
| Analysis documentation | R Markdown |

## Limitations

- **Validation leakage:** Differential-expression filtering uses all samples before the initial train/test split. The resulting LASSO feature list is reused in later cross-validation and repeated hold-outs instead of being reselected within each training fold. The reported AUCs therefore remain exploratory; these experiments do not establish generalization to an independent cohort.
- **Metric interpretation:** `Control` is the first factor level, so `twoClassSummary` sensitivity refers to controls ([caret documentation](https://topepo.github.io/caret/measuring-performance.html)). The custom `QLIKE` formula differs from binary log loss. The confusion-matrix plot replaces computed counts with a [hard-coded matrix](final%20coding.Rmd#L1199-L1227), so its counts are not used as verified results here.
- **Stability comparison:** The [hold-out versus final-model plot](final%20coding.Rmd#L1159-L1185) repeats each final AUC 20 times; these are not 20 independent final-model evaluations.
- **Feature rankings:** The notebook averages absolute LASSO coefficients with `caret` importance scores, which are on different scales. SVM/kNN importance from `caret` uses a model-independent filter rather than the permutation procedure described in the notebook ([caret documentation](https://topepo.github.io/caret/variable-importance.html)).

## Repository Structure

```text
.
├── README.md          # Project summary
└── final coding.Rmd   # Complete analysis source
```

## How to Run

**Preparation is required; a complete HTML render has not been verified.**

1. Install the dependencies below.
2. Obtain `GSE68801_series_matrix.txt.gz` from the linked GEO record, then replace the absolute Windows path in `getGEO(filename = ...)` with your local path.
3. Open `final coding.Rmd` in RStudio and execute chunks in order. Check the metadata first: the initial code uses `gender:ch1 == "M"`, while the repeated hold-out block uses `Sex:ch1 == "male"`. These fields and values must be reconciled against the actual data before a full run.

<details>
<summary>Dependency setup</summary>

```r
install.packages(c(
  "BiocManager", "R.utils", "reshape2", "ggplot2", "dplyr", "tidyr",
  "pheatmap", "RColorBrewer", "DT", "caret", "randomForest",
  "glmnet", "pROC", "e1071", "kernlab", "rmarkdown"
))
BiocManager::install(c("GEOquery", "limma", "hgu133plus2.db"))
```

This includes packages loaded by the original notebook and the SVM engine dependency. Package versions are not pinned. Another source inconsistency: the earlier tuning note states kNN `k = 9`, while the final model code uses `k = 3`.

</details>

## Key Takeaways

The recorded analysis suggests a disease-associated signal in selected expression features, with LASSO showing the strongest reported final classification result. Its value as a portfolio project is the connected workflow: preparing metadata, working with a high-dimensional matrix, combining statistical filtering with regularization, comparing classifiers, and communicating results through diagnostic plots and feature interpretation. The validation limitations also show why judging the boundaries of the evidence matters.
