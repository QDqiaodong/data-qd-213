# data-qd-213 · 光谱波长校准

利用已识别的参考谱线拟合二次像素到波长关系，保留逐线残差、整体验收状态，并只在标定范围内换算下游像素。

这是实验数据处理 API，不做峰检测、元素识别或图形展示。字段、精度和失败边界见 docs/contract.md。

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
