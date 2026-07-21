# 投稿与返修可重复性核对表

## 首次投稿前

- 在最终代码状态运行 `Rscript scripts/run_all.R --force`，或在已校验原始数据上运行完整离线重建。
- 确认 `results/reproducibility/validation_report.csv` 每一行 `passed = TRUE`。
- 保存 `artifact_inventory.csv`、`last_pipeline_run.csv`、`environment/package_versions.csv`、`environment/system_info.csv` 和 `environment/sessionInfo.txt`。
- 确认 Figure 1–5 及 Figure S5 均有同名 SVG、PDF、600-dpi TIFF、PNG preview 和独立 source-data CSV。
- 打开 PDF/SVG 检查字体、特殊字符、面板标记、轴标题和置信区间；不要只检查 PNG 缩略图。
- 在方法学和图注中明确 baseline OR、RCT interaction、纵向变化差，以及 sample-first、donor-inference、Match=Yes donor-site 配对规则。
- 核对 Figure 5 只含定位、donor map 和关键基因；基线缓解、聚合敏感性与纵向结果只出现在 Figure S5。
- 核对 Figure S5 明确报告 legacy 上皮方向与 bulk 不一致（OR 2.90，FDR 0.123），且 sample-first 结果未被用于重新选基因。
- 把 `config/data_manifest.csv` 与 `config/predefined_gene_sets.csv` 随代码归档，确保输入和主签名可追溯。

## 返修修改后

- 只修改审稿意见涉及的脚本或配置；不要手工修改 `data/derived/` 或 `results/`。
- 用 `Rscript scripts/run_all.R --from=<受影响的最早阶段>` 重建所有下游产物。
- 对比新旧 `artifact_inventory.csv`，确认变化只出现在预期文件。
- 对比 `table_cross_mechanism_gate.csv`、`table_longitudinal_gate.csv` 与单细胞三项 gate，若任何预设判定改变，先停止写作并重新评估主线。
- 更新图注中的 n、效应值、95% CI 和探索性标签；不复制旧版数字。
- 将返修运行的 `logs/pipeline/`、环境快照和 validation report 独立归档。

## 可接受的重建差异

同一原始文件、冻结基因集和包环境下，CSV 数值与图形应确定性重建。不同操作系统的字体替代可能造成文字换行或极小的栅格像素差异，但不应改变表格数值、面板数据或科学判定门。若 MD5、样本数、供者数、效应方向或预设 gate 不一致，应视为失败而不是“正常波动”。
