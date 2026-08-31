# `bm_pirna` 脚本说明

## `find_homolog_pairs.py`

`find_homolog_pairs.py` 使用 DIAMOND 对两个蛋白 FASTA 文件执行双向
`blastp`，并在基因层面筛选同源基因对。脚本默认采用双向最佳命中
（reciprocal best hit，RBH）模式，默认筛选阈值如下：

- 氨基酸一致性（identity）≥ 40%；
- 查询序列和目标序列的覆盖度均 ≥ 70%；
- E-value ≤ 1e-10。

蛋白到基因的映射可以自动从 NCBI FASTA 标题中的 `[GeneID=...]`、
FlyBase FASTA 标题中的 `parent=FBgn...`，以及蛋白下载流程生成的
`protein_download_manifest.tsv` 中解析。也可以通过 `--map1` 和 `--map2`
提供显式的蛋白—基因映射文件。脚本支持多异构体：多个蛋白可以映射到
同一基因，并为每个基因对保留最佳蛋白比对结果。

`--mode` 支持以下配对模式：

- `rbh`（默认）：保留基因层面的双向最佳命中；
- `best`：保留第一个 FASTA 中每个基因在第二个 FASTA 中的单向最佳命中；
- `all`：保留所有通过筛选阈值的非重复基因对。

输出 TSV 以 `gene1`、`gene2` 和 `homology`（identity 百分比）为核心列，
同时包含 `gene1_label`、`gene2_label`、`protein1`、`protein2`、
`alignment_length`、`query_coverage_pct`、`subject_coverage_pct`、`evalue`、
`bitscore` 和 `pair_method`。其中 `protein1` 和 `protein2` 记录产生该基因对
最佳比对结果的蛋白 accession。

覆盖度参数支持 `--query-cover`/`--min-query-cover` 和
`--subject-cover`/`--min-subject-cover` 两组别名。每次运行还会在结果文件
旁生成 `<输出文件名>.metadata.json`，记录输入文件校验和、DIAMOND 版本、
运行参数、映射来源及命中统计。
