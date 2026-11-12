# 超阈值尾部分析方法与输出

`python -m climate_data.cli tail <资料.json> <配置.json>` 对逐日雨量做超阈值（POT）尾部分析，输出单个JSON报告。资料口径见[资料约定](data-contract.md)；本文件说明方法、求解与输出字段。配置示例见`configs/tail-normal.json`。

## 分析前给定的参数

- `u`：阈值（mm，非负），严格超过才计为超阈值日。
- `run_length`：runs法去簇参数（非负整数）。
- `min_segment_days`：最短完整段长度（公历日数）。
- `sigma_bounds`/`xi_bounds`：GPD拟合的公开有限矩形，sigma为正。
- `return_periods`：重现期T列表（均>1）。
- `sensitivity`：一组`{u, run_length}`，与主配置共用其余参数，全部重算并保留。

## 段审核与暴露

声明观测段内缺行或null即整段拒绝（不补干日），短于`min_segment_days`的段同样拒绝；重叠段、重复日期、负雨量在读取层拒绝。暴露=所有合格完整段日数/365.2425，包含段边界观察日。被拒绝的段连同原因（`incomplete`/`below-min-segment-days`及缺日计数）写入报告。

## runs法去簇与边界簇

严格超过u的日按runs法成簇：相邻超阈值日之间至少有`run_length`个连续已观测非超阈值日才分成两簇。簇峰取簇内最大值，同峰取最早日。簇前后若不足`run_length`个连续已观测非超阈值日（按段内位置计），按边界簇排除并写明原因；已隔开的内部簇不因边缘簇被排除而失效。每个簇（入选与排除）都报告首尾日、超阈值日数、峰日、峰值、超量及前后非超阈值连续日数，事件独立性近似可逐簇追查。边界排除只是估计约定，不声称消除全部相关性。

## GPD拟合与求解诊断

对入选簇峰超量z拟合GPD：NLL在xi趋零时用稳定级数（指数极限），支持域`1+xi*z/sigma>0`为硬约束，越界返回inf而非容差忽略。求解为矩形约束Nelder-Mead多起点（15个矩形分位点+1个PWM点），另做41×41粗网格扫描。报告包含：

- `starts`：每个起点的起止点、NLL、收敛标记与函数评估次数；
- `grid_scan`：可行网格点数与网格最优；
- `distinct_local_optima`：去重后的局部最优及命中次数；
- `candidate`：全部起点中NLL最小者；
- `kkt`：候选点解析梯度、参数边界状态（interior/lower-bound/upper-bound）、KKT符号检查、支持余量`support_min_margin`与所用容差。

`solver_status`恒为`local-candidates-only`：全局最优未证明，候选只是矩形内多起点与网格扫描的最佳局部结果，这一局限写入`limitations`。

## 重现水平

lambda=入选簇数/暴露年数。按Poisson簇到达+GPD尾部，年最大值不超概率为`exp(-lambda*S(z))`，令其为`1-1/T`得所需生存概率`S=-log(1-1/T)/lambda`：

- `S>1`：目标水平在阈值u以下，`status=below-threshold`，超出尾部适用域，不给水平值、不外推；
- 否则由逆函数求z与水平`u+z`；xi<0时报告有限端点与端点余量，余量非正记`endpoint-violation`；
- 每个T都独立做往返检查：由水平回算生存概率与年不超概率，报告残差（`roundtrip`）。

## 报告结构

顶层含`source_sha256`（输入身份）、`config_sha256`、量纲`units`、原始摘要`raw_summary`、回显参数`parameters`、公开约定`conventions`、局限`limitations`、`primary`与`sensitivity`（每组含`u`/`run_length`回显、段审核、暴露分母、入选簇、排除边界簇、GPD拟合、lambda与重现水平）。无入选簇或无可行候选时对应字段为null并附`note`。敏感性各组结果全部保留，未择优剔除。
