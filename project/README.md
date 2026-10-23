# YZDHD AD disease-side analysis

035最终科学复核已于2026-10-11完成，覆盖15靶点、29个背景组合，状态为`FINAL_SCIENTIFIC_REVIEW_COMPLETE_WITH_LIMITATIONS`。数值审计通过，未新增已确认机制靶点；细胞身份、低表达匹配和治疗方向的限制保留，药物等级仍为C0。见[最终报告](results/2026-10-11-virtual-ko-scientific-review-035/最终科学复核报告.md)、[逐组合结论](results/2026-10-11-virtual-ko-scientific-review-035/target_verification_conclusions.tsv)和[复现入口](src/run_final_scientific_review_035.py)。

本项目沿用既有 G1–G5 冻结定义。G2c 人类 H4K12la 转译补证按 2026-09-30 修订 016 收束，G6 空间验证按 2026-09-29 修订 015 收束。现阶段目标是整理 G1 人类 AD 转录关联与 G2–G5 模型范围证据，形成供未来益智地黄汤真实扰动数据比较的有方向状态接口。

## 当前锁定范围

- 研究数据资源：D1–D5、D7–D8、G5 的 GSE261698、GSE291239、GSE74614、GSE290101；G6 数据集仅保留历史筛查审计记录，D4 仅用于描述性空间背景。
- G1 主发现/验证队列：GSE160936、GSE148822、GSE174367；D6/GSE214979 保留为探索性队列。
- G5 主队列：GSE261698；验证与敏感性队列不能挽救主门控失败。
- G2c：人类转译补证已关闭，当前工作状态为 `G2C_CLOSED_INSUFFICIENT_HUMAN_EVIDENCE`；历史结果 `HUMAN_TRANSLATION_UNRESOLVED` 与总体 `G2_PASS_MODEL_SCOPED` 保留。
- G6：已关闭，状态为 `G6_CLOSED_NO_QUALIFIED_DIRECT_ACCESS_DATA`；历史执行状态 `G6_NOT_EXECUTABLE` 保留。筛过的数据集没有通过全部 G6 门槛；没有计算 AD–对照空间效应。逐路线处置见 `results/2026-09-29-g6-closure-015/`。
- 主比较：AD vs control；连续病理指标作为敏感性分析。
- G1–G5 统计单位沿用各自冻结协议；当前项目不进行 G6 空间推断。
- 数据策略：优先使用公共 processed 数据，原始数据保持只读。
- 药物侧：按修订 017 部分收束，状态为 `DRUG_VALIDATION_PARTIALLY_CLOSED_MISSING_DATA`，证据上限 C0。缺失数据依赖的模型新验证、药物状态恢复、暴露—效力、整方药效/机制及药物介导通讯验证退出当前待执行范围；保留先验、G1–G5 来源登记、疾病接口整理与现有证据可完成的治疗方向审计。未知方向/参照保持待定。MAP 已运行但药物输出不变，排序不放行，修复不在当前范围。当前入口见 [药物侧当前状态](docs/protocol/drug_side_latest.md)，逐项处置及重启要求见 [收束报告 017](results/2026-09-30-drug-validation-partial-closure-017/report.md)。
- 固定处方：修订018登记用户提供的16味药名；名称级公开数据库成分整理、结构/ADMET计算初筛已在020/021完成，剂量、炮制、制备和批次信息不是该生信分析的前置条件。022 校正 ChEMBL 相似结构误列为直接靶点记录的问题；023 澄清与 G1 稳健244基因的交集：导入的数据库靶点映射先验有8个探索性交集，但只覆盖6/2,276个候选CID；1条待查原始assay的结构匹配BindingDB记录交集为0。方剂实测成分/暴露未建立，药物等级仍为C0，017关闭的实验验证项未重启。见[交集口径复核023](results/2026-09-30-fixed-formula-g1-overlap-scope-023/report.md)。

固定处方公共靶点补查024已执行：[正式SOP](docs/protocol/2026-09-30-fixed-formula-target-supplement-SOP-024.md)、[最新结果与覆盖](results/2026-09-30-fixed-formula-target-supplement-024/report.md)、[机器审计](results/2026-09-30-fixed-formula-target-supplement-024/machine_audit.json)。新实验记录、辅助实验、预测、相似结构及待核验关联独立报告；017关闭范围、C0及G1–G5结论保持。022/023统计为历史输入范围，不替代024。

## 目录

- `data/raw/`：原始下载文件与不可修改的来源记录。
- `data/interim/`：格式转换、基因映射、QC 中间文件。
- `data/processed/`：供下游分析使用的冻结输入。
- `src/`：正式分析代码。
- `results/`：表格与图形输出。
- `docs/protocol/`：方案与实施协议。
- `docs/decisions/`：关键研究决策记录。
- `params.yaml`：影响结果的集中参数配置。

## 复现入口

项目保留 Level 2 的 Git、参数配置、Decision Log 和环境锁定要求。G6 输入分支已退役，目前没有活跃 G6 数据处理 stage；raw 数据保持只读，DVC remote 尚未验证。

## 当前状态

G6 历史记录仍为 `G6_NOT_EXECUTABLE`，未计算疾病组间效应；修订 015 已结束当前空间验证路线，原因是没有候选通过全部门槛，不代表空间效应阴性。G5 已完成并判定为 `G5_PASS_MODEL_SCOPED`，结果见 `results/2026-09-22-disease-side-v3.4-g5/`；该结论限于已测试的 APP/PS1-21 动物星形胶质细胞响应和模型内证据。G2c 历史结果仍为 `HUMAN_TRANSLATION_UNRESOLVED`；修订 016 已因现有证据不足关闭该补证路线，并未将其改判为生物学阴性。项目不据此声称已在人脑验证空间关联、细胞间通信或 H4K12la 因果机制。

GSE233208 的历史 DVC 输入预检已完成，随后该路线按既有决策退役；当前没有活跃 G6 stage。该预检不运行 RCTD、不计算组间效应。G2c 与 G6 的当前状态、关闭依据和重启条件见 `docs/protocol/disease_side_latest.md`。
