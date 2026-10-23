param([string]$ProjectRoot = (Split-Path $PSScriptRoot -Parent))
$ErrorActionPreference = 'Stop'
$base = Join-Path $ProjectRoot 'results/2026-10-04-stage-b-tcmsp-docking-032'
$source = Join-Path $base 'audited_docking_results.tsv'
$auditPath = Join-Path $base 'output_audit.json'
$rows = @(Import-Csv -LiteralPath $source -Delimiter "`t")
$audit = Get-Content -LiteralPath $auditPath -Raw | ConvertFrom-Json
$manifest = @(Import-Csv (Join-Path $base 'ligand_manifest.tsv') -Delimiter "`t")
if ($rows.Count -ne 23 -or ($rows.ligand_id | Sort-Object -Unique).Count -ne 23) { throw 'Task count or uniqueness mismatch' }
if (@(Compare-Object ($rows.ligand_id | Sort-Object) ($manifest.ligand_id | Sort-Object)).Count) { throw 'Manifest mismatch' }
if (($rows | Measure-Object n_poses -Sum).Sum -ne 230 -or $audit.poses -ne 230) { throw 'Pose count mismatch' }
$labels = @{'969516'='姜黄素'; '1794427'='绿原酸'; '164676'='丹参酮IIA'; '5281794'='6-姜烯酚'}
$detail = foreach ($r in $rows) {
    $label = if ($r.role -eq 'positive_control') { if ($r.target -eq 'AKR1B1') {'tolrestat（对照）'} else {'ARM1（对照）'} } elseif ($labels.ContainsKey($r.cid)) { $labels[$r.cid] } else { "CID$($r.cid)" }
    $note = if ($r.role -eq 'positive_control') {'对照姿态回收有支持；不抵消受体准备疑点'} elseif ($r.target -eq 'LTA4H') {'条件候选；Zn保留及来源真实性未解决；不能推断双酶选择性'} elseif ($r.cid -eq '164676') {'沿用原审计：优先在保留NADP条件下复核；未确认抑制效力'} elseif ($r.cid -in @('969516','1794427')) {'沿用原审计：最佳姿态与原始NADP严重重叠，暂不接受排名'} else {'保留待复核；结合辅因子距离及逐原子几何检查，未作新阈值筛选'}
    [pscustomobject][ordered]@{target=$r.target; label=$label; cid=$r.cid; ligand_id=$r.ligand_id; role=$r.role; database_herb_association_unverified=$r.herbs; source_chemical_name=$r.names; input_state=$r.input_state; input_charge=$r.input_charge; output_charge=$r.output_charge; score_kcal_mol=$r.score_kcal_mol; poses=$r.n_poses; min_cofactor_distance_A=$r.min_distance_to_requested_cofactor_A; exact_input_identity_preserved=$r.exact_input_identity_preserved; interpretation=$note}
}
$detail | Export-Csv -LiteralPath (Join-Path $base '对接任务汇总.csv') -NoTypeInformation -Encoding utf8BOM
$candidates = @($detail | Where-Object role -eq 'candidate')
$groups = @($candidates | Group-Object target,cid)
if ($groups.Count -ne 16 -or $candidates.Count -ne 20) { throw 'Candidate counts mismatch' }
$compound = foreach ($g in $groups) {
    $items = @($g.Group | Sort-Object {[double]$_.score_kcal_mol})
    $first = $items[0]
    [pscustomobject][ordered]@{target=$first.target; label=$first.label; cid=$first.cid; task_count=$items.Count; best_score_kcal_mol=$first.score_kcal_mol; all_state_scores=($items | ForEach-Object { "$($_.input_state):$($_.score_kcal_mol)" }) -join '; '; min_cofactor_distances_A=($items.min_cofactor_distance_A -join '; '); database_herb_association_unverified=$first.database_herb_association_unverified; source_chemical_name=$first.source_chemical_name; interpretation=$first.interpretation}
}
$compound = @($compound | Sort-Object target,{[double]$_.best_score_kcal_mol})
$compound | Export-Csv -LiteralPath (Join-Path $base '候选成分汇总.csv') -NoTypeInformation -Encoding utf8BOM
$lines = [System.Collections.Generic.List[string]]::new()
$lines.Add('# 当前分子对接结果汇总')
$lines.Add('')
$lines.Add('整理日期：2026-10-04。范围为已查收的032批次；仅重整已有结果与审计，不重新对接、不新增评分阈值、不改变药物侧C0证据等级。')
$lines.Add('')
$lines.Add('**23/23任务完成，230个构象；涉及16种不同候选成分和2种阳性对照。AKR1B1为15种候选成分、19份候选状态输入及2份对照输入；LTA4H为1种候选成分及1份对照输入。输入状态不能当作独立成分或独立重复。计算完成，正式结构验收未通过。**')
$lines.Add('')
$lines.Add('## 候选成分完整表')
$lines.Add('')
$lines.Add('按靶点分别展示，表内排序仅是最佳Vina评分顺序，不是药效或复核优先级。单位：评分kcal/mol；距离Å。未确认的常用名称保留CID，原始系统名称完整保留在CSV。药材关联均来自数据库，不能认定为煎液实测来源。')
$lines.Add('')
$lines.Add('|靶点|成分/CID|各输入状态评分|与原始辅因子最近重原子距离|解释|')
$lines.Add('|---|---|---|---|---|')
foreach ($r in $compound) { $lines.Add("|$($r.target)|$($r.label) / $($r.cid)|$($r.all_state_scores)|$($r.min_cofactor_distances_A)|$($r.interpretation)|") }
$lines.Add('')
$lines.Add('## 阳性对照')
$lines.Add('')
$lines.Add('|靶点|对照|输入状态|评分|最佳分数姿态共晶RMSD（Å）|')
$lines.Add('|---|---|---|---:|---:|')
foreach ($r in @($detail | Where-Object role -eq 'positive_control')) {
    $control = @($audit.controls | Where-Object ligand_id -eq $r.ligand_id)
    if ($control.Count -ne 1) { throw 'Missing control RMSD' }
    $rmsd = ([double]$control[0].best_score_pose_rmsd_A).ToString('F3',[cultureinfo]::InvariantCulture)
    $lines.Add("|$($r.target)|$($r.label)|$($r.input_state)|$($r.score_kcal_mol)|$rmsd|")
}
$lines.Add('')
$lines.Add('## 关键审计与收敛建议')
$lines.Add('')
$lines.Add('1. AKR1B1导出复合物21/21无NAP，LTA4H导出复合物2/2无Zn；实际prepared receptor未交付，尚不能区分计算时遗漏和仅导出时遗漏。最高分AKR候选与原始NADP的极端重叠否定其最佳姿态在该辅因子存在时的几何可行性。')
$lines.Add('2. AKR1B1—丹参酮IIA优先复核；姜黄素、绿原酸暂不接受最高分排名，但保留在修正受体后重新搜索。其余候选保留，不能仅凭距离没有极端重叠认定合格。LTA4H—6-姜烯酚仍为条件候选，方剂来源和双酶功能问题未解决。')
$lines.Add('3. 5份输入电荷/精确身份被平台改变；中性与羧酸根输入部分折叠到同一输出状态，状态间分差不能解释为pH效应。')
$lines.Add('4. 对照RMSD支持姿态回收，但不能证明整个受体准备流程有效。版本、搜索seed、exhaustiveness、energy_range与实际受体文件缺失，不能严格重现云端计算。4Å接触清单只是邻域描述，不等于氢键。')
$lines.Add('5. 当前首先应核验实际受体与辅因子，再复核少量组合；本次未把疾病关联、数据库药材来源或对接分数升级为药效、直接靶点或表观遗传因果证据。')
$lines.Add('')
$lines.Add('## 文件、方法与复现')
$lines.Add('')
$lines.Add('- 23任务完整表：对接任务汇总.csv；16候选成分合并表：候选成分汇总.csv。')
$lines.Add('- 上游：audited_docking_results.tsv、output_audit.json、ligand_manifest.tsv、结果解读.md；原始结构和230姿态沿用既有存档。')
$lines.Add('- 受体：AKR1B1/2FZD、LTA4H/4L2L；24×24×24Å盒，中心与准备参数见configs/stage_b_tcmsp_docking_032.json。准备seed42不是云端搜索seed。')
$lines.Add('- 本次在本机PowerShell用原生CSV/JSON操作汇总；复用已有审计结果和实验技术库22_docking_simulation证据边界，无新增依赖、远程计算或模型训练。')
$lines.Add('- 复现：在项目根目录执行 pwsh -File src/summarize_stage_b_docking_032.ps1。23任务唯一性、输入清单一致性、230构象、16候选成分及3份对照映射均在脚本中检查；不重算原审计RMSD。')
$lines.Add('- 本次为可追溯的确定性结果整理；上游仅达到探索性验证，尚未达到严格计算复现。DVC remote恢复未在本次验证。无一次性脚本/缓存产生；本脚本为正式复现资产。')
$lines | Set-Content -LiteralPath (Join-Path $base '当前对接结果汇总.md') -Encoding utf8
$provenance = [ordered]@{date='2026-10-04'; mode='deterministic aggregation of existing audit'; task_count=$rows.Count; candidate_tasks=$candidates.Count; distinct_candidate_compounds=$groups.Count; poses=230; source_sha256=(Get-FileHash -LiteralPath $source -Algorithm SHA256).Hash; audit_sha256=(Get-FileHash -LiteralPath $auditPath -Algorithm SHA256).Hash; script_sha256=(Get-FileHash -LiteralPath $PSCommandPath -Algorithm SHA256).Hash; validation='task uniqueness, manifest equality, pose counts, compound counts, control mappings passed'; new_docking=$false; new_model_training=$false}
$provenance | ConvertTo-Json | Set-Content (Join-Path $base 'summary_provenance.json') -Encoding utf8
Write-Output 'PASS: 23 tasks, 20 candidate tasks, 16 compounds, 230 poses; summary and two CSVs generated.'
