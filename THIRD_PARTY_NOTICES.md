# 第三方材料与引用

## OCx24 数据

Jehad Abed, Jiheon Kim, Muhammed Shuaibi, Brook Wander, Boris Duijf, Suhas Mahesh, Hyeonseok Lee, Vahe Gharakhanyan, Sjoerd Hoogland, Erdem Irtem, Janice Lan, Niels Schouten, Anagha Usha Vijayakumar, Jason Hattrick-Simpers, John R. Kitchin, Zachary W. Ulissi, Aaike van Vugt, Edward H. Sargent, David Sinton, C. Lawrence Zitnick. **Open Catalyst Experiments 2024 (OCx24): Bridging Experiments and Computational Models.** 2024. [arXiv:2411.11783](https://arxiv.org/abs/2411.11783), [DOI](https://doi.org/10.48550/arXiv.2411.11783).

[官方数据页面](https://facebookresearch.github.io/fairchem/ocx24/) 明确将数据集按 [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) 提供。本发布包以 FAIR-Chem 提交 **0a625a8924cf0a87dffbcd2f0ac58d82ca43bfa4** 为唯一计算来源，不跟随上游主分支自动更新。

随包上游数据：ExpDataDump_241113.csv、ExpDataDump_241113_clean.csv、XRFDataDump-241113.csv、XRDDataDump-241113.csv、supporting_data/batches_info.csv、processed_data/CO2R_40_70_all.csv。其路径、固定版本原始下载 URL、字节数及 SHA256 见 [下载清单](01_sources/download_manifest.json)。这些 CSV 字节保持原样；`02_data/` 中的配对、解析、分组和质量标记属于本项目派生处理，不是 OCx24 作者新提供的数据。

本项目的推断与评价不代表 OCx24 作者、Meta 或其机构的认可。报告中的固定版本单位审计也不宣称已核实当前上游版本的状态。

## FAIR-Chem 上游代码

仅附一份用于阅读和单位审计的 `ocx/core/data_handling.py`，不执行或导入该代码。其固定版本路径为 [FAIR-Chem / applications / ocx](https://github.com/facebookresearch/fairchem/tree/0a625a8924cf0a87dffbcd2f0ac58d82ca43bfa4/src/fairchem/applications/ocx)。

Copyright (c) Meta, Inc. and affiliates. MIT License，完整文本保留于 [01_sources/LICENSE.md](01_sources/LICENSE.md)。项目原创代码在根目录 MIT 许可下单独署名。

## 建议同时引用

1. 本二次分析：JIN HANZHE. *A grouped evaluation of HER premeasurement for CO selectivity prediction in OCx24*. v1.0.0，本地发布候选版，2026。机器可读信息见 [CITATION.cff](CITATION.cff)。正式发布后应补充实际仓库与版本链接。
2. 上述 Abed 等的 OCx24 论文及固定版本数据来源。

本包没有附其他论文全文或第三方付费附件。

