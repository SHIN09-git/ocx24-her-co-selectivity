# OCx24：HER 前测对 CO 选择性预测的增量检验

[English](README.en.md) · [完整中文报告](05_reports/首期计算研究报告.md) · [复现说明](REPRODUCING.md) · [署名](AUTHORSHIP.md) · [许可](LICENSING.md)

**在固定的 OCx24 数据、输入和模型范围内，加入 HER 前测没有带来稳定的 CO 选择性预测收益。** 这是一项可复现的负结果与方法研究，不是新材料发现，也不能推出 HER 在所有条件下都没有价值。

作者：**JIN HANZHE**，Independent Researcher, Chongqing, China。联系：jinhanzhe09@gmail.com。原分析日期：2026-10-12；公开版整理日期：2026-10-22。版本：**[v1.0.0](https://github.com/SHIN09-git/ocx24-her-co-selectivity/releases/tag/v1.0.0)（2026-10-22）**。未进行独立人工科学复核或同行评审。

## 研究问题与主要结果

在 CO₂ 测试开始前，测得电极在 50 mA/cm² 下的 HER 全电池电压，能否改善其在 200 mA/cm² 下 CO 法拉第效率（CO FE）的预测？法拉第效率表示电量用于生成某种产物的比例。

M0 使用目标配方的元素比例与合成来源；M1 仅增加同电极 HER 电压的绝对值。主算法为受限随机森林，另用岭回归与 XRF 组成子集检查稳健性。误差按制备编号等权计算，单位为**百分点**，不是百分比相对变化。

|主比较：配方 + 随机森林|电极 / 制备编号 / 元素组合|M0 MAE|M1 MAE|改善 M0−M1|改善的 95% 区间|
|---|---:|---:|---:|---:|---:|
|开发集：元素组合留出|404 / 216 / 38|13.76|13.80|−0.04|[−1.79, 1.44]|
|原始预留评价|37 / 14 / 10|11.26|12.33|−1.06|[−3.16, 0.28]|

表内数值分别四舍五入；改善使用未舍入的误差计算。区间来自固定预测下的元素组合重采样，不包含所有训练与实验不确定性。

公开价值在于：保留负结果；提供从原始行号到配对、分组、预测、评分的完整证据链；展示重复电极、批次差异、组成外推和筛选命中率怎样影响评价。方法可作为相似小样本材料数据的审计例子。

## 两种复现入口

从 [v1.0.0 发布页](https://github.com/SHIN09-git/ocx24-her-co-selectivity/releases/tag/v1.0.0) 下载源代码包并解压，或克隆本仓库。先安装 Python 3.12，再进入包含复现入口的目录。**快速核查只用 Python 标准库，不需要下载数据或安装第三方包：**

```console
python reproduce.py
```

它检查发布文件哈希，核对 441 个原始有效配对、20 条指定记录、6 项映射异常、216 个批次日期差异，并独立计算 23 项已保存结果的统计。它不训练模型，也不等于重新得到全部预测。

**从随包原始 CSV 开始完整复算：**

```console
python -m pip install -r requirements.txt
python reproduce.py --full
```

完整复算重新构建数据、拟合原定模型、执行对照和原始预留顺序，逐项比较保存结果。两种模式均在新的 `reproductions/` 子目录输出日志和机器可读报告，不覆盖发布证据。也可用 `--output 新目录` 指定位置；已有目录会被拒绝。仅检查文件完整性可用 `python reproduce.py --verify-only`。

依赖安装需要网络；安装完成后的复算使用随包数据，可离线运行。详细步骤、输出范围与实测验证见 [复现说明](REPRODUCING.md) 和 [发布验证记录](VALIDATION.md)。

## 阅读入口

- [完整中文报告](05_reports/首期计算研究报告.md) / [浏览器阅读版](05_reports/report.html)：方法、结果、图表与局限。
- [数据字典](05_reports/数据字典.md)：字段含义、单位、分组和质量标记。
- [冻结方案](00_admin/preregistered_analysis_plan.json)、[原始划分](00_admin/split_manifest.json)、[原始解封记录](00_admin/holdout_unsealed.json)：保留历史字节与时间。
- [开发集结果](04_results/primary_comparisons.csv)、[预留结果](04_results/holdout_results.json)：主要数值。
- [上游来源清单](01_sources/download_manifest.json)、[打包变更记录](provenance/PACKAGING.md)：版本、下载来源与哈希。
- [引用信息](CITATION.cff)、[发布说明](RELEASE_NOTES.md)：供后续公开使用。

## 结论的边界

共 441 个电极配对、230 个制备编号、48 个元素组合；制备编号不自动等同独立合成。预留集是在公开数据预审计之后划分的，原方案只在本地冻结。文件名中的 `preregistered` 不表示已公开预注册。本公开版中的标签已知，复算是历史分析的重放，不能称为新的盲测或独立验证。

HER 前测不是随机处理，现有数据不能回答“先做 HER 是否会改变后续 CO₂ 性能”。来源留出也同时改变元素覆盖，不能把它当作实验室总体的独立外部验证。当前结果不支持继续原计划的正向材料筛选路线；未执行扩展模型或湿实验。

报告的乙烯摩尔产率附录只审计下述固定上游版本；它不表示乙烯 FE 错六倍，不影响本项目的 CO FE 主结论，也不宣称上游当前版本仍有同一问题。

## 数据来源、署名与许可

原始数据来自 [Abed 等的 OCx24](https://arxiv.org/abs/2411.11783)，固定于 FAIR-Chem 提交 `0a625a8924cf0a87dffbcd2f0ac58d82ca43bfa4`。随包只保留复算所需的六份数据表、一份上游代码和上游许可文件，未改写这些原始文件。数据集的许可见 [OCx24 官方说明](https://facebookresearch.github.io/fairchem/ocx24/)。

本项目原创代码采用 MIT；原创报告、图表及派生表采用 CC BY 4.0。上游数据保留其 CC BY 4.0 署名要求，上游代码保留 Meta 的 MIT 许可。使用时应同时引用本项目与 OCx24；具体范围见 [LICENSING.md](LICENSING.md) 和 [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md)。

来源核查、计算、文字和公开版整理使用了 OpenAI Codex；AI 不列为作者。详见 [AUTHORSHIP.md](AUTHORSHIP.md)。

