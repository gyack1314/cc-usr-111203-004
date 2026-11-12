# 暴雨超阈值尾部的暴露量与重现水平

这是独立的气候资料Python项目。`data/`保存合成资料，`climate_data.io`负责严格读取、公历日期和来源字节摘要，`climate_data.domain`提供现有描述统计，`climate_data.cli`将摘要输出为JSON。所有观测值、站点、人群、情景与金额都是虚构。字段和研究约定见[资料约定](docs/data-contract.md)。

当前Python 3.11及以上可直接使用，没有第三方Python依赖。项目根目录运行：

```sh
python -m unittest discover -s tests -v
python -m compileall -q climate_data tests
python -m climate_data.cli inspect data/normal.json
```

`inspect`显示当前原始摘要，不覆盖资料。`invalid-`前缀文件为不合法资料，其他变体使用相同字段保留业务边界。来源字节SHA-256随摘要输出，后续结果应同时标明输入身份、量纲和适用范围。

`tail`运行超阈值尾部分析：runs法去簇（含边界簇排除）、GPD多起点矩形约束拟合与KKT诊断、Poisson簇重现水平及往返检查，并对一组(u,run_length)做敏感性重算，全部结果写入单个JSON报告：

```sh
python -m climate_data.cli tail data/normal.json configs/tail-normal.json
```

分析配置（阈值、run_length、最短完整段长度、GPD矩形边界、重现期、敏感性组）在分析前给定，示例见`configs/tail-normal.json`；方法与输出字段见[尾部分析说明](docs/tail-analysis.md)。

容器使用固定基础镜像，可运行`docker build -t climate-material .`，随后`docker run --rm climate-material`检查资料。`docker run --rm climate-material python -m unittest discover -s tests -v`执行已有检查。运行不依赖远端账号或实时服务。
