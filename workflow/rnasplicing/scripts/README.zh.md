# RNA 剪接流程脚本说明

本目录包含 RNA 剪接分析流程中使用的 3 个辅助脚本。

## 概览

- **`group_files.py`**
  - 根据 `config.json` 和 STAR 比对结果生成样本分组文件与比较索引文件。

- **`summarize_splicing.py`**
  - 读取单个 comparison 的 rMATS 结果，合并不同事件类型，导出汇总表，并绘制火山图。

- **`analyze_clusters.py`**
  - 基于 `featureCounts` 输出分析 piRNA cluster 表达变化，计算 RPK、log2FC，并生成汇总图。

---

## 1. `group_files.py`

### 功能

扫描 STAR 生成的 BAM 文件目录，根据 `config.json` 中定义的样本分组生成下游 rMATS 所需的分组文件和比较索引文件。

### 主要输入

- **`--bam-dir`**
  - STAR 比对结果目录。
  - 目录结构示例：

```text
star/
  CONTROL_REP1/
    Aligned.sortedByCoord.out.bam
  CONTROL_REP2/
    Aligned.sortedByCoord.out.bam
  SUGP1_REP1/
    Aligned.sortedByCoord.out.bam
```

- **`--config`**
  - JSON 格式的流程配置文件。
  - 该脚本主要使用以下字段：

```json
{
  "samples": {
    "control": ["CONTROL_REP1", "CONTROL_REP2"],
    "sugp1": ["SUGP1_REP1", "SUGP1_REP2"]
  },
  "comparisons": [
    ["control", "sugp1"]
  ]
}
```

- **`--output-dir`**
  - 输出目录。

### 主要输出

- **`{output_dir}/{group}.txt`**
  - 每个样本组一个文件。
  - 格式：单行、多个 BAM 路径以逗号分隔。

```text
star/CONTROL_REP1/Aligned.sortedByCoord.out.bam,star/CONTROL_REP2/Aligned.sortedByCoord.out.bam
```

- **`{output_dir}/comparison_{i}_{group1}_vs_{group2}.txt`**
  - 每个比较一个文件。
  - 格式：一行，两列，使用制表符分隔。

```text
control	sugp1
```

- **`{output_dir}/comparison_index.txt`**
  - 所有比较的索引文件。
  - 格式：每行一个 comparison。

```text
comparison_1_control_vs_sugp1	control	sugp1
```

### 示例

```bash
python group_files.py \
  --bam-dir star \
  --config config.json \
  --output-dir groups
```

---

## 2. `summarize_splicing.py`

### 功能

读取单个 comparison 下的 rMATS 结果表，合并支持的剪接事件类型，输出汇总表，并基于 `IncLevelDifference` 和 `FDR` 绘制火山图。

该脚本还支持读取带有 `is_known` 和 `ENTREZ.ID` 列的注释 CSV，对已知 piRNA 相关基因进行打标签。

### 支持读取的 rMATS 事件文件

脚本会在 `--input-dir` 下查找以下文件：

- `A3SS.MATS.JC.txt`
- `A5SS.MATS.JC.txt`
- `MXE.MATS.JC.txt`
- `RI.MATS.JC.txt`
- `SE.MATS.JC.txt`

### 主要输入

- **`--input-dir`**
  - 单个 comparison 的 rMATS 结果目录。

- **`--output-file`**
  - 合并汇总表输出路径。

- **`--plot-file`**
  - 火山图 PDF 输出路径。

- **`--fdr`**
  - 点着色使用的 FDR 阈值。
  - 默认值：`0.05`

- **`--annotation-file`**
  - 可选的注释 CSV 文件，用于打标签。
  - 预期列包括：

```text
sample_id,silkbase_id,ENTREZ.ID,Gene_Symbol,is_known,Gene_Name
```

### 输入表格式要求

每个 rMATS 结果表至少应包含这些列：

```text
ID
GeneID
geneSymbol
PValue
FDR
IncLevel1
IncLevel2
IncLevelDifference
```

其中 `GeneID` 可以是：

- 数值型 ID，例如 `101741153`
- 带 `LOC` 前缀的 ID，例如 `LOC101741153`

如果提供了 `--annotation-file`，脚本会从 `LOC...` 中提取数值部分，再与 `ENTREZ.ID` 进行匹配，只对 `is_known=true` 的行打标签。

### 主要输出

- **`summary.txt` / 自定义汇总路径**
  - 制表符分隔的合并结果表。
  - 输出列为：

```text
ID	GeneID	geneSymbol	PValue	FDR	IncLevel1	IncLevel2	IncLevelDifference
```

  - 字段含义：
    - `ID`
      - rMATS 事件 ID。
      - 在单个事件文件内通常唯一，用于标识一条剪接事件记录。
    - `GeneID`
      - 与该剪接事件关联的基因 ID。
      - 可能是纯数值 ID，也可能是 `LOC` 前缀形式，例如 `LOC101741153`。
    - `geneSymbol`
      - 基因符号。
      - 如果 rMATS 原始结果中缺失，该列可能为空。
    - `PValue`
      - 该剪接事件差异检验的原始 P 值。
    - `FDR`
      - 多重检验校正后的显著性值。
      - 火山图纵轴使用的是 `-log10(FDR)`。
    - `IncLevel1`
      - 第 1 组样本的 inclusion level 值。
      - 通常为逗号分隔的多个重复值。
    - `IncLevel2`
      - 第 2 组样本的 inclusion level 值。
      - 通常为逗号分隔的多个重复值。
    - `IncLevelDifference`
      - 两组 inclusion level 的差值，通常可理解为 `group2 - group1` 的 ΔΨ。
      - 火山图横轴使用该列。

- **`volcano.pdf` / 自定义 PDF 路径**
  - 火山图，包含：
    - 灰色：不显著事件
    - 蓝色：显著下降事件
    - 红色：显著上升事件
    - 可选：known piRNA 基因标签

### 示例

```bash
python summarize_splicing.py \
  --input-dir rmats/comparison_1_control_vs_sugp1 \
  --output-file results/comparison_1_control_vs_sugp1/summary.txt \
  --plot-file results/comparison_1_control_vs_sugp1/volcano.pdf \
  --annotation-file ../../data/external/piRNA_annotations_id.csv
```

---

## 3. `analyze_clusters.py`

### 功能

基于 `featureCounts` 输出的 count matrix 分析 piRNA cluster 表达变化。脚本会进行长度归一化得到 RPK，计算分组均值 log2 fold change，评估重复间一致性，并输出统计表与图形结果。

### 主要输入

- **`--counts`**
  - `featureCounts` 输出表。

- **`--config`**
  - 包含样本分组定义的 `config.json`。

- **`--group1`**
  - 参考组名称。

- **`--group2`**
  - 处理组名称。

- **`--output-dir`**
  - 结果输出目录。

- **`--min-rpk`**
  - 用于过滤的最小平均 RPK。
  - 默认值：`1.0`

- **`--log2fc`**
  - 图中使用的绝对 log2FC 阈值。
  - 默认值：`1.0`

- **`--pseudocount`**
  - log2 变换前添加的伪计数。
  - 默认值：`0.5`

### 输入表格式要求

脚本期望输入标准 `featureCounts` 结果，包含：

- 元信息列：
  - `Geneid`
  - `Chr`
  - `Start`
  - `End`
  - `Strand`
  - `Length`
- 每个 BAM 一个 count 列

表头示例：

```text
Geneid	Chr	Start	End	Strand	Length	star/CONTROL_REP1/Aligned.sortedByCoord.out.bam	star/CONTROL_REP2/Aligned.sortedByCoord.out.bam
```

### 主要输出

- **`{group2}_vs_{group1}_clusters.tsv`**
  - 制表符分隔结果表。
  - 包含：
    - cluster ID
    - 两组的平均 RPK
    - 平均 log2FC
    - 每个重复的 log2FC 列，如 `log2FC_REP1`、`log2FC_REP2`
  - 常见字段含义：
    - `cluster`
      - piRNA cluster 的 ID，对应输入表中的 `Geneid`。
    - `mean_rpk_{group1}`
      - 参考组所有重复的平均 RPK。
    - `mean_rpk_{group2}`
      - 处理组所有重复的平均 RPK。
    - `mean_log2FC`
      - 两组平均表达的 log2 fold change，计算方式为 `log2((group2 + pseudocount) / (group1 + pseudocount))`。
    - `log2FC_REP1`, `log2FC_REP2`, ...
      - 按重复配对后得到的每对样本 log2 fold change，用于评估重复一致性。

- **`{group2}_vs_{group1}_rep_correlation.pdf`**
  - 重复间 log2FC 相关性散点图。
  - 仅在至少存在两对重复时生成。

- **`{group2}_vs_{group1}_log2fc_distribution.pdf`**
  - 平均 log2FC 分布直方图。

### 示例

```bash
python analyze_clusters.py \
  --counts results/cluster_counts.txt \
  --config config.json \
  --group1 control \
  --group2 sugp1 \
  --output-dir results/cluster_analysis
```

---

## 说明

- `group_files.py` 和 `summarize_splicing.py` 使用 `argparse`。
- `analyze_clusters.py` 使用 `typer`。
- 图形输出均为 PDF 文件。
- `summarize_splicing.py` 当前只读取 `JC` 模式的 rMATS 文件。
