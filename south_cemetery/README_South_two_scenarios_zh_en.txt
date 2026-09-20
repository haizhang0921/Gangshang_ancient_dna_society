Gangshang South Cemetery chronology — two kinship scenarios
2026-09-20

一、分析目的
-----------
本工作包用于重新计算岗上南区墓地的年代模型，并对 SM1–SM2 亲属关系的不确定性进行敏感性分析。

所有情景使用完全相同的 OxCal archaeology-only 输出 gangshang_s.csv。该 OxCal 模型已经包含：
1) SM14 与 SM9 之间的考古先后关系，通过 Sequence() 表示；
2) SM1 中四个个体作为同一顶层 burial/depositional event，通过 Combine() 表示；
3) Cemetery Start / Cemetery End 以及 Cemetery Duration 的 OxCal 输出。

R 第二阶段不会再次乘入 SM14–SM9 的顺序约束，也不会把 SM1R1–SM1R4 当作独立墓葬重新输入，因此避免重复使用同一考古信息。

二、两种亲属关系情景
-----------------
Main（正文）:
  genetic_relations_s.csv
  SM1–SM2 = 2nd relationship
  absolute death/burial-event gap = 35 ± 26 years
  direction established = No

Alternative（附录敏感性分析）:
  genetic_relations_s_alt.csv
  SM1–SM2 = Sibling
  absolute death/burial-event gap = 26 ± 22 years
  direction established = No

由于南区只有一条 kinship edge，SM1 与 SM2 的二维联合后验可以在 OxCal 五年网格上精确枚举，不需要 MCMC。其余 12 个顶层 burial-event marginal posteriors 在 kinship update 中保持不变。

三、顶层年代事件
-------------
程序应识别 14 个顶层 burial events：
SM13, SM10, SM14, SM9, SM17, SM15, SM3, SM12, SM6, SM2, SM1, SM8, SM11, SM5。

SM1R1–SM1R4 是 Combine("SM1") 的内部测年成分，不重新进入第二阶段模型。
总计应识别 17 个 R_Date determinations + 1 个顶层 Combine event SM1。

四、运行方法
-----------
将 R working directory 设置为本文件夹，然后：

  source("run_two_scenarios.R")

程序依次运行：
1) Main exact kinship model
2) Main conditional boundary/duration reconstruction
3) Alternative sibling exact kinship model
4) Alternative conditional boundary/duration reconstruction
5) 两情景比较

主要结果目录：
  south_R_exact_main/
  south_R_boundaries_main/
  south_R_exact_alt_sibling/
  south_R_boundaries_alt_sibling/

主要比较文件：
  south_two_scenarios_burial_mean_comparison.csv
  south_two_scenarios_relationship_diagnostics.csv
  south_two_scenarios_boundary_duration_comparison.csv
  south_two_scenarios_key_results.csv

五、Boundary reconstruction 的性质
--------------------------------
与北区一样，OxCal CSV 只提供 marginal posterior distributions，并不提供原始 OxCal full joint MCMC states。因此，R 中对 Start / End / Duration 的更新属于 conditional modular reconstruction，而不是对原始 OxCal Boundary posterior 的精确恢复。

边界重建使用 14 个 burial-event joint draws。SM1–SM2 由精确二维 posterior 联合抽样，其余 burial events 从各自 OxCal marginals 抽样。条件 Phase likelihood 使用 uniform event-intensity 形式。程序内部沿用经过验证的 cal BP 计算形式，同时额外输出 *_calBCE.csv 供论文使用；duration 始终以 years 表示。

六、出版约定
-----------
- Conventional radiocarbon determinations: BP ± laboratory error.
- Calibrated/modelled calendar ages: cal BCE/CE.
- Durations: years.
- 不把 posterior ranges 称为 confidence intervals；使用 posterior probability / HPD terminology。

七、独立核对值
-------------
reference_exact_pair_expected.csv 是根据上传的 gangshang_s.csv 对 SM1–SM2 二维 posterior 做独立精确计算得到的核对值，不是正式 R 输出。正式论文结果应以 R 输出为准。

对当前输入，预计：
- Main 2nd-degree 情景：SM1 posterior mean ≈ 2527.24 cal BCE；SM2 ≈ 2533.45 cal BCE；mean absolute gap ≈ 34.59 years。
- Alternative sibling 情景：SM1 posterior mean ≈ 2526.94 cal BCE；SM2 ≈ 2530.90 cal BCE；mean absolute gap ≈ 28.65 years。

这说明 sibling alternative 预计主要改变 SM2 和 SM1–SM2 gap，对整个南区墓地时间框架的影响应很小；最终以正式 R boundary/duration 输出核定。

------------------------------------------------------------
English analytical note
------------------------------------------------------------
The two scenarios use the same archaeology-only OxCal marginals. The SM14–SM9 Sequence and the SM1 Combine event are already represented in those marginals and are not imposed a second time in R. The only difference between scenarios is the undirected SM1–SM2 temporal-gap likelihood: 35 ± 26 years for the main second-degree model versus 26 ± 22 years for the sibling sensitivity model. Because the South Cemetery contains only one kinship edge, the SM1–SM2 joint update is evaluated exactly on the finite OxCal grid rather than by MCMC.
