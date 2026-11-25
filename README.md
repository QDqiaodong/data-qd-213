# data-qd-213 · 光谱波长校准

利用已识别的参考谱线拟合二次像素到波长关系，保留逐线残差、整体验收状态，并只在标定范围内换算下游像素。

这是实验数据处理 API，不做峰检测、元素识别或图形展示。字段、精度和失败边界见 docs/contract.md。

## 计数守恒的产品更正与可追溯交付

`wavelength_calibration.SpectrumArchive(root)` 在原 `process` 精确拟合之上提供持久档案：校准血缘按整数版本更正（旧版证据永久可读），原始曝光登记后不可变，`build` 在固定波长网格上精确分摊计数（目标箱计数之和 + dropped_counts 恒等于原始计数），`release` 以 request_id 幂等交付当前 ready 产品并冻结证据快照。契约细节见 docs/archive_contract.md，可复现示例见 examples/archive_sample.py。

## 安装与验证

Python 3.11+，运行只用标准库；构建固定 setuptools 80.9.0、wheel 0.45.1。

    uv venv --python 3.12 .venv
    uv pip install --python .venv/bin/python --offline -r requirements-build.lock
    uv pip install --python .venv/bin/python --offline --no-index --no-build-isolation -e .
    .venv/bin/python -m unittest discover -s tests -v
    .venv/bin/python examples/run_sample.py

样例曲线为 500 + pixel + pixel²/1000，参考线残差为零，查询 50/150/250 像素得到可复核波长。

## 边界

本版要求调用方已经完成参考线识别，且像素和波长同向递增。不处理多谱阶、波长反转、异常点剔除或浮点协方差。
