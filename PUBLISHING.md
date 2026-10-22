# 发布与维护记录

公开仓库：[SHIN09-git/ocx24-her-co-selectivity](https://github.com/SHIN09-git/ocx24-her-co-selectivity)  
正式版本：[v1.0.0](https://github.com/SHIN09-git/ocx24-her-co-selectivity/releases/tag/v1.0.0)  
发布日期：2026-10-22。

版本附件提供完整源代码包、可单独阅读的中文 HTML 报告及 SHA256 校验和。仓库中的 `RELEASE_MANIFEST.json` 覆盖发布文件，清单本身除外。

## 后续维护顺序

1. 在新提交中修改目标文件，保留原冻结证据。若修改计算代码、输入或模型配置，重新运行完整复算。
2. 如修改中文报告，安装 `requirements-report.txt` 后运行 `python tools/build_report.py` 更新 HTML。
3. 更新版本、引用文件和发布说明，再运行 `python tools/make_manifest.py` 冻结新清单。该工具只用于维护者有意修改后的版本，不用于掩盖损坏检查失败。
4. 运行快速核查和适用的完整复算，检查最终源代码包能否解压并复现。核对线上工作流对应的提交；状态入口见 [Actions](https://github.com/SHIN09-git/ocx24-her-co-selectivity/actions/workflows/reproduce.yml)。
5. 将确定的提交标记为新版本，附带源代码包、报告和新校验和。已发布的 v1.0.0 标签和附件作为历史版本保留。

只修改网址或发布文字，也应重新检查引用、文档链接、文件哈希与压缩包。保留原始完整重训记录，公开标签不会重新成为盲测数据。

