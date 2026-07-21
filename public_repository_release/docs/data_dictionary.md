# 数据字典

本字典覆盖从队列整理到 Figure 1–5 的主要分析对象。机器可读版本为 `docs/data_dictionary.csv`。原始 GEO/H5AD 文件中的完整发布者字段保持原样；这里只定义进入本项目估计量、判定门和图形的数据。

## 通用编码与统计方向

- `response_binary = 1` 表示达到队列预先定义的临床应答或缓解终点，`0` 表示未达到。
- `remission = Remission` 为 GSE282122 未来缓解；`Non_Remission` 为未缓解。
- 所有 `*_z` 评分均在每个 accession 的基线 UC 分布内标准化。不同平台的原始 GSVA 数值不直接混合。
- `Inflammation_z` 为冻结的 `HALLMARK_INFLAMMATORY_RESPONSE` 主评分；基因列表见 `config/predefined_gene_sets.csv`。
- `delta_* = follow-up − baseline`。炎症评分为负表示分子状态下降/逆转。
- 基线模型的 `odds_ratio` 是评分每升高 1 SD 的应答 OR；小于 1 表示高炎症与较低应答相关。
- 纵向 `estimate` 是 responder 变化减去 nonresponder 变化；小于 0 表示 responder 下降更深。
- 单细胞 `module_score` 使用 baseline Inflamed、sample-first、供者内样本等权的 donor-state 参考集对三个区室共同基因中心化和标准化，再对有效冻结基因取均值。
- 置信区间均为 95%；缺失区间通常表示供者数或结局事件数未达到预设最低要求，而不是效应为零。

## 核心数据层

| 数据层 | 文件 | 行单位 | 主要用途 |
|---|---|---|---|
| 原始输入清单 | `config/data_manifest.csv` | 一个公开文件 | 下载、来源和 MD5 校验 |
| 队列注册表 | `data/derived/extended/sample_registry.csv` | 一个 GEO 样本 | 队列整理与纳排追踪 |
| 基线患者数据 | `results/source_data/common_state_baseline_patient_data.csv` | 一个 subject×cohort×arm | 基线 Firth 模型与 RCT 分臂模型 |
| 协变量增强基线数据 | `results/source_data/covariate_augmented_baseline_patient_data.csv` | 一个 subject×cohort | 仅公开且可无歧义匹配的临床协变量 |
| 协变量可用性矩阵 | `results/tables/common_state/table_covariate_availability_matrix.csv` | cohort×requested covariate | 公开范围、缺失和可调整性审计 |
| 有限调整模型 | `results/tables/common_state/table_adjusted_baseline_associations.csv` | 一个治疗队列 | 未调整与调整后炎症评分 OR 并列 |
| 机制间差异检验 | `results/tables/common_state/table_mechanism_moderator_test.csv` | 一个分析集 | 完整与调整子集的 REML/Hartung–Knapp omnibus 检验与机制 OR 幅度 |
| 分层异质性 | `results/tables/common_state/table_stratified_sensitivity_meta.csv` | 一个分层 | 终点、时间、机制和研究层面活检部位敏感性 |
| 样本排除审计 | `results/tables/common_state/table_sample_exclusion_audit.csv` | cohort×互斥原因 | 从原始 accession 样本到独立参与者的流程 |
| 纵向配对数据 | `results/source_data/longitudinal_paired_patient_data.csv` | 一个 subject×cohort 配对 | 分子变化与高状态逆转 |
| 单细胞样本-状态 pseudobulk | `results/source_data/singlecell_sample_state_pseudobulk.csv` | sample×state×compartment | 单细胞第一聚合层 |
| 单细胞样本-区室 pseudobulk | `results/source_data/singlecell_sample_compartment_pseudobulk.csv` | sample×compartment | 临床审计第一聚合层 |
| 单细胞定位 pseudobulk | `results/source_data/singlecell_compartment_pseudobulk.csv` | donor×state×compartment | baseline Inflamed、样本等权的 Figure 5 推断输入 |
| 单细胞精确纵向配对 | `results/source_data/singlecell_longitudinal_site_pairs.csv` | donor×site×compartment | Match=Yes 的 Post−Pre 配对 |
| 单细胞定位表 | `table_singlecell_compartment_localization.csv` | compartment×state | donor 均值与置信区间 |
| 图源数据 | `results/source_data/Figure1–5_source_data.csv` | panel 对应的最小绘图行 | 投稿 source-data 交付 |

## 单细胞聚合规则

三个 h5ad 分别是 `epicolonic_final.h5ad`、`myeloid_final.h5ad` 和 `fibperi_final.h5ad`。仅纳入 UC、`Pre`/`Post`、结局为 `Remission`/`Non_Remission` 且元数据完整的细胞。表达矩阵按 CSR 非零值分块读取，先形成至少 20 个细胞的 sample×state 与 sample×compartment pseudobulk；同一供者的多个活检样本等权。主定位限定 baseline Inflamed 活检。纵向仅纳入 `Match=Yes`，按 donor×site 精确匹配 Pre/Post，再在供者内对部位等权。旧 donor×state×time/细胞数加权结果仅作为明确标记的 legacy sensitivity 保留。完整规则见日期锁定的单细胞分析修正说明。

## 关键结果字段

详细逐字段定义见 CSV。投稿表述时尤其注意：`arm` 区分 Active/Placebo；`estimand` 明确当前行是 OR、交互还是变化差；`fdr` 只在对应分析族内校正；`baseline_high` 是队列内基线最高三分位，不是临床切点；`molecular_reversal` 是随访评分低于该队列基线中位数。
