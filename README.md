<p align="center">
  <img src="https://raw.githubusercontent.com/21century-Aissa/Sirius/main/SIRIUS_dialogs/sirius_logo.png" alt="SIRIUS Logo" width="180"/>
</p>

# AI Sirius MOFA — Intelligent Agent for Multi-Omics Analysis and Precision Medicine

AI Sirius MOFA is an intelligent bioinformatics agent designed to automate and assist multi-omics analysis in the context of precision medicine. It combines a desktop graphical interface (CustomTkinter) with an R-based analytical pipeline built on MOFA2 (Multi-Omics Factor Analysis v2), a Bayesian latent factor model that simultaneously integrates multiple omics layers — mutations, methylation, mRNA, miRNA, proteomics — to identify coordinated and interpretable sources of biological variation. The agent is powered by an LLM (Ollama, Claude, or DeepSeek) that guides the user at every step of the analysis, interprets results, diagnoses pipeline errors, and automatically searches for evidence in the PubMed literature.

### AI Sirius MOFA — 面向多组学分析与精准医学的智能代理

AI Sirius MOFA 是一款智能生物信息学代理，专为精准医学场景下的多组学分析自动化与辅助而设计。它将桌面图形界面（CustomTkinter）与基于 MOFA2（多组学因子分析 v2）的 R 分析流程相结合。MOFA2 是一种贝叶斯潜因子模型，可同时整合多个组学层——突变、甲基化、mRNA、miRNA、蛋白质组——以识别具有生物学意义的协同变异来源。该代理由大语言模型（Ollama、Claude 或 DeepSeek）驱动，在分析的每个步骤为用户提供指导，解释结果，诊断流程错误，并自动检索 PubMed 文献证据。

---

## Installation

Before launching AI Sirius MOFA, ensure the following dependencies are installed on your system.

On the R side, install the required packages from an R console:
```r
install.packages(c("survival", "survminer", "msigdbr"))
BiocManager::install("MOFA2")
BiocManager::install("ReactomePA")
```

Make sure `Rscript` is accessible in the system PATH (verify with `Rscript --version` in a terminal). On the Python side, install the project dependencies (from the package root folder):
```bash
pip install -r requirements.txt
```

If you are using Ollama as a local LLM provider, start the daemon before launching the application:
```bash
ollama serve
```

For Claude or DeepSeek, set the corresponding environment variable (`ANTHROPIC_API_KEY` or `DEEPSEEK_API_KEY`) before starting the application. Once all prerequisites are met, launch the interface:
```bash
python ai_sirius_mofa.py
```

The health bar in the interface automatically checks the availability of Rscript and the LLM at startup.

### 安装

在启动 AI Sirius MOFA 之前，请确保系统中已安装以下依赖项。

在 R 环境中，从 R 控制台安装所需包：
```r
install.packages(c("survival", "survminer", "msigdbr"))
BiocManager::install("MOFA2")
BiocManager::install("ReactomePA")
```

确保 `Rscript` 已添加到系统 PATH 中（在终端运行 `Rscript --version` 验证）。在 Python 环境中，从项目根目录安装依赖：
```bash
pip install -r requirements.txt
```

如果使用 Ollama 作为本地 LLM 提供商，请在启动应用前先启动守护进程：
```bash
ollama serve
```

如果使用 Claude 或 DeepSeek，请在启动前设置对应的环境变量（`ANTHROPIC_API_KEY` 或 `DEEPSEEK_API_KEY`）。完成所有前置条件后，启动界面：
```bash
python ai_sirius_mofa.py
```

界面的健康状态栏会在启动时自动检测 Rscript 和 LLM 的可用性。

---

## Data Preparation

Each omics matrix must follow this convention: rows correspond to features (genes, CpG sites, proteins), columns to samples (patient identifiers identical across all views), and values must be numeric. Recommended formats are: `log2(TPM+1)` for mRNA and miRNA, beta values `[0,1]` for methylation, binary `0/1` matrix for mutations (bernoulli likelihood), and `log2 LFQ` for proteomics. The optional metadata file must have its first column named exactly `sample`, with identifiers matching the omics matrix column names exactly (case-sensitive). Additional clinical variables (condition, survival time, event, subtype) are added as extra columns.

### 数据准备

每个组学矩阵必须遵循以下规范：行对应特征（基因、CpG 位点、蛋白质），列对应样本（所有视图中患者标识符必须完全一致），值必须为数值型。推荐格式为：mRNA 和 miRNA 使用 `log2(TPM+1)`，甲基化使用 beta 值 `[0,1]`，突变使用二进制 `0/1` 矩阵（伯努利似然），蛋白质组使用 `log2 LFQ`。可选的元数据文件的第一列必须命名为 `sample`，其中的标识符须与组学矩阵的列名完全匹配（区分大小写）。其他临床变量（状态、生存时间、事件、亚型）作为额外列添加。

---

## Usage — MOFA2 Pipeline

From the graphical interface, select your omics files in the corresponding tabs, define the working directory (`work_dir`), choose the number of factors (recommended starting value: 15), and launch the main pipeline. The `mofa_pipeline.R` script runs in the background via `Rscript`; R progress and logs are visible in real time in the Log panel of the interface. The trained model is saved as `<work_dir>/model.hdf5` and can be reloaded at any time with `loadModel()`.

Once the model is trained, the Visualization, Survival, and GSEA tabs allow you to generate additional plots and analyses. Available visualization types include: data overview (`data_overview`), variance explained per factor and per view, feature weight heatmaps, factor scatter plots, and biplots. The Survival tab runs Cox analyses (summary, hazard ratios as forest plots) and Kaplan-Meier curves stratified on the most significant factor. The GSEA tab performs gene set enrichment analysis on selected factors via ReactomePA and msigdbr.

### 使用 — MOFA2 流程

在图形界面中，在对应标签页选择您的组学文件，定义工作目录（`work_dir`），选择因子数量（推荐初始值：15），然后启动主流程。`mofa_pipeline.R` 脚本通过 `Rscript` 在后台运行；R 的进度和日志在界面的 Log 面板中实时显示。训练完成的模型保存为 `<work_dir>/model.hdf5`，可随时通过 `loadModel()` 重新加载。

模型训练完成后，可通过可视化、生存分析和 GSEA 标签页生成更多图表和分析。可用的可视化类型包括：数据概览（`data_overview`）、各因子和各视图的解释方差、特征权重热图、因子散点图和双标图。生存分析标签页支持 Cox 分析（汇总、森林图展示风险比）以及基于最显著因子分层的 Kaplan-Meier 曲线。GSEA 标签页通过 ReactomePA 和 msigdbr 对所选因子进行基因集富集分析。

---

## Usage — AI Assistant

The conversational assistant is permanently accessible in the Chat panel. It operates in three modes automatically detected from your question: `MOFA_ASSISTANCE` for interpreting results and pipeline usage guidance, `PIPELINE_ERROR_HELP` for diagnosing R errors (it analyzes the error text, identifies the likely cause, and proposes the minimal fix), and `PUBMED_EVIDENCE` for searching literature evidence — enable the PubMed option in settings for the agent to automatically query NCBI E-utilities and cite PMIDs in its responses.

The assistant automatically loads text files present in the `work_dir` (CSV, TSV, TXT, LOG, PDF) to ground its answers in your actual results rather than generic knowledge. It maintains a short conversation history (last 12 exchanges) to ensure coherent dialogue without overloading the LLM context. You can switch LLM provider at any time from the settings; the conversation history is automatically reset upon switching.

### 使用 — AI 助手

对话助手可随时在 Chat 面板中访问。它根据您的问题自动检测并切换三种模式：`MOFA_ASSISTANCE` 用于解释结果和流程使用指导，`PIPELINE_ERROR_HELP` 用于诊断 R 错误（分析错误文本，识别可能原因并提出最小化修复方案），`PUBMED_EVIDENCE` 用于检索文献证据——在设置中启用 PubMed 选项后，代理将自动查询 NCBI E-utilities 并在回答中引用 PMID。

助手会自动加载 `work_dir` 中的文本文件（CSV、TSV、TXT、LOG、PDF），以便将回答建立在您的实际结果之上，而非泛泛的通用知识。它维护一个简短的对话历史（最近 12 条交互），在不超载 LLM 上下文的前提下保持对话连贯性。您可以随时在设置中切换 LLM 提供商，切换时对话历史将自动重置。
