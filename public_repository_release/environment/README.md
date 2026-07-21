# 运行环境说明

分析入口不会自动安装或升级任何包。首次配置运行 `Rscript scripts/install_dependencies.R`，随后用 `Rscript scripts/00_preflight.R` 核验并冻结本机实际环境快照。

## 参考环境

- R 4.5.1（本次开发与核查环境）
- macOS / Unix-like shell；其他平台应提供等价的 `bash` 与 `curl`
- CRAN：`logistf`, `metafor`, `lme4`, `ggplot2`, `patchwork`, `ggridges`, `ggrepel`, `scales`, `svglite`, `ragg`, `dplyr`
- Bioconductor：`GEOquery`, `Biobase`, `AnnotationDbi`, `hthgu133pluspm.db`, `hgu133plus2.db`, `GSVA`, `rhdf5`
- R 推荐内存 16 GB；完整原始数据和中间产物建议预留至少 20 GB 磁盘

`Matrix` 随 R 库安装并由 preflight 明确核验。`msigdbr` 只用于基因集清单缺失时的人工引导式重建；正式流水线读取版本化的 `config/predefined_gene_sets.csv`，因此主运行不依赖在线基因集版本。

## 每次运行生成的环境证据

`scripts/00_preflight.R` 写出：

- `package_versions.csv`：每个必需/可选包的可用状态和精确版本；
- `system_info.csv`：R、平台、操作系统、时区、locale、项目路径和命令行工具路径；
- `sessionInfo.txt`：完整 R session 信息。

各分析脚本还在 `logs/` 写入独立 `sessionInfo`。返修时应连同 `config/data_manifest.csv`、`config/predefined_gene_sets.csv`、`results/reproducibility/validation_report.csv` 一起归档。

## 大型 h5ad 的资源控制

单细胞脚本按 CSR 非零值分块读取，默认每块 5,000,000 个非零值；同一遍读取同时构建 sample×state、sample×compartment 和 legacy 审计聚合。内存受限时可缩小：

```sh
H5AD_CHUNK_NNZ=1000000 Rscript scripts/run_all.R --from=singlecell_analysis
```

改变块大小不会改变数值定义。每个区室聚合完成后会写入带输入 MD5、冻结签名 MD5 和 sample-first 算法版本的 RDS 缓存；只有三者完全一致时才复用。

## 跨平台注意事项

正式图、预览和视觉核查均由 R 导出：SVG 使用 `svglite`，PDF 使用 base R vector device（`useDingbats=FALSE`），TIFF/PNG 使用 `ragg`。字体默认 Helvetica。本机 R 的 Cairo/XQuartz 动态库不可用，因此不调用 `cairo_pdf`；这不改变 R-only 后端或可编辑矢量文本。Linux 无 Helvetica 时会使用系统映射字体，几何和统计结果不变，但提交前应检查 PDF 字体嵌入和换行。MD5 在 R 内用 `tools::md5sum()` 核验；单细胞下载脚本会自动选择 macOS 的 `md5` 或 Linux 的 `md5sum`。
