# 计数守恒的光谱产品更正与可追溯交付契约

本文档定义 `wavelength_calibration.SpectrumArchive`。它在既有 `process` 精确二次拟合之上，提供校准血缘版本化、原始曝光登记、固定波长网格上的计数守恒产品，以及面向下游的幂等交付。`process` 本身的输入输出与成功、拒绝、异常行为不变。

运行只依赖 Python 3.11+ 标准库；交付为库 API，不含命令行工具、页面或绘图。只处理已识别参考线和单调递增的单谱段：不做峰检测、谱线识别或异常点剔除。

## 数值与身份契约

- 所有外部身份（`slot_id`、`exposure_id`、`request_id`、`release_id`）为非空字符串。
- 新 API 的数值输入是长度不超过 32 的普通十进制字符串（可选符号、整数或小数）。拒绝 NaN、无穷、分数字面量（`"1/2"`）、科学计数（`"1e3"`）以及 Python 数字（`int`/`float`）。
- 内部一律用精确有理数计算，不存在“浮点近似后再转分数”。数值输出统一为最简 `n` 或 `n/d` 字符串；返回值为可 JSON 序列化的普通 dict/list/str/int。
- 返回对象与归档内部状态隔离：调用方修改任何返回层级的列表或字典，都不会改变归档。

## 档案

`SpectrumArchive(root, *, fault_hook=None)`，`root` 是调用方提供的独立持久化目录。首次打开空目录创建档案；关闭（`close()`）后重新打开同一目录可继续工作。不依赖进程全局变量、工作目录、网络或宿主已有数据。同一目录不要求支持多进程同时操作，也不抵御磁盘硬件损坏。所有成功返回的写入在返回前已持久化（原子替换 + fsync）。

## 校准血缘：`set_calibration` / `get_calibration`

`set_calibration(slot_id, expected_revision, data)`：

- `data` 含旧契约的 `reference_lines`（3..1000 条、像素与波长同向严格递增）与 `max_residual_nm`。调用方不提供 `query_pixels`；实现以参考线两端像素作为查询调用原 `process`，把系数、逐线残差、最大残差、验收状态原样存为证据。拟合能力完全来自现有 `process`，没有第二套拟合口径。
- 创建血缘时 `expected_revision=0`；更正时必须等于该血缘当前整数版本。成功后版本加 1（从 1 开始）。版本不符、输入非法、残差 rejected 均抛 `ValueError`：不保存候选、不消耗版本、不影响已可用产品。
- 除旧契约外，新 API 额外要求：拟合多项式在整个参考像素闭区间上导数严格大于 0，且映射波长为正。导数是线性函数，检查两端即可判定整个闭区间；不满足即拒绝。例如参考像素 0/1/2 对应波长 400/401/404（容差 0）旧拟合 accepted，但拟合为 `400 + p²`，在像素 0 处导数为 0，必须拒绝。旧 `process` 的行为不因此改变。
- 每次成功更正都是新证据版本：即使只改参考线 id、系数残差与最终计数完全不变，也不与旧版合并。
- 返回记录含 `slot_id`、`revision`、原始输入 `input`（`reference_lines` 与 `max_residual_nm`）、`coefficients`、`reference_residuals`、`max_abs_residual_nm`、`status`、`calibrated`（两端点查询证据）。

`get_calibration(slot_id, revision=None)`：省略 `revision` 取当前版，指定版本取当时记录；查询不重拟合、不覆盖历史。未知 slot 或版本抛 `ValueError`。

## 原始曝光：`add_exposure`

`add_exposure(exposure_id, slot_id, pixel_edges, counts, wavelength_edges)` 登记不可变原始曝光并返回所登记记录。

- slot 必须已存在；`exposure_id` 不可复用，即使内容相同也抛 `ValueError`。
- 设 N 个源箱：`pixel_edges` 长 N+1，`counts` 长 N；M 个目标箱：`wavelength_edges` 长 M+1。1<=N,M<=2000。边界严格递增；像素边界非负、波长边界为正、计数非负。
- 登记不检查当前校准是否覆盖边界；覆盖资格在 `build` 时判断，以便保留暂时无法处理的原始曝光。登记后原始曝光、目标网格与所属 slot 不可修改。

## 产品：`build` / `inspect_exposure`

`build(exposure_id)` 用该曝光所属 slot 的当前版本生成产品：

- 所有源像素边界必须落在参考像素范围内；所有目标波长边界必须落在映射后的源波长范围内。不外推、不把未观测部分补零；任一条件不满足抛 `ValueError`，不留下半个新产品。目标范围可以严格小于源范围。
- 算法：每个源像素箱两端由二次多项式映射为波长区间，该箱积分计数在映射区间内均匀分布；目标箱获得每个源箱计数乘以“波长交叠长度 / 该源箱完整波长宽度”，再精确求和。不是把像素中心映射后投到最近目标箱，也不是按像素长度分摊。边界左闭右开、最后右端闭合；零长度接触不贡献计数。
- 产品含 `product_id`、`exposure_id`、`slot_id`、`calibration_revision`、原始像素边界与计数、目标波长边界、映射后的源波长边界 `source_wavelength_edges`、逐目标箱计数 `target_counts`、目标范围外的 `dropped_counts`，以及完整校准输入与拟合证据 `calibration`。
- 精确满足 `sum(target_counts) + dropped_counts == sum(counts)`。`dropped_counts` 是目标范围之外的源计数，不是把覆盖不足隐藏成丢弃——覆盖不足在 build 时直接失败。
- 同一曝光在同一当前校准版已有成功产品时，重复 `build` 原样返回旧产品，不新增身份；成功更正后即使数值相同，下一次 `build` 也必须产生新的 `product_id` 与证据版本。旧产品不可改写。

`inspect_exposure(exposure_id)` 返回原始曝光及 `state`、`last_product`：

- `unbuilt`：从未成功 build；`ready`：最后成功产品的 `calibration_revision` 等于当前版；`stale`：存在最后成功产品但版本不符。`last_product` 在没有时为 `null`，否则为完整旧产品。
- 成功的校准更正立即使该 slot 全部已生成曝光变为 `stale` 并撤销其新交付资格；其他 slot 不受影响。资格判定只看血缘版本匹配，不用“重新计算后计数未变”代替；重建一个曝光不会连带恢复同 slot 其他曝光。build 失败保留最后产品与 stale 状态。

## 交付：`release` / `get_release`

`release(request_id, exposure_ids)` 把一批不同曝光的当前 ready 产品交付给下游：

- 输入必须非空、无重复、无未知曝光，且任一曝光 unbuilt/stale 时整批失败：不交付子集、不静默 rebuild、不占用 `request_id`、不消耗号码。失败请求随后补齐产品后，可用同一 `request_id` 正常成功。
- 成功时返回含唯一 `release_id`、连续整数 `release_number`（从 1 开始，只对成功的新请求加 1，失败、读取、重放不产生缺号）、`request_id` 和按 `exposure_id` 稳定排序的完整产品快照列表。
- 成功交付冻结每份产品的原始数据、网格、校准身份与逐线证据；之后的更正、重建或新交付都不改变它。`get_release(release_id)` 始终返回当时完整交付——它不代表这些旧产品现在仍获准用于新交付。
- `request_id` 是跨关闭/重开的幂等身份，请求内容为曝光 id 的无序集合（不含调用当时的当前版本）。合法且内容相同的已成功请求优先原样重放旧交付（旧 `release_id`、`release_number`、`product_id` 与证据），即使曝光后来 stale 或已 rebuild，也不改交付、不重新授权、不再计号。同一 `request_id` 用于不同曝光集合抛 `ValueError`。非空身份、非空列表、重复项检查先于重放判断，重复项不会被无序集合折叠。只有第一次成功的新请求才按当前资格取产品。

## 中断边界与 fault_hook

校准当前版本切换与所有相关曝光资格变化一起生效；整批 release 的产品快照、`request_id` 绑定与计号一起生效。重新打开档案不会看到“已变更当前版但旧曝光仍 ready”“有交付但重试又生成第二份”“占了请求或号码但缺交付”的中间状态。此义务同样适用于 `add_exposure` 与 `build` 的身份和完整数据。

`fault_hook(phase, operation)` 用于测试这些边界：`operation ∈ {set_calibration, add_exposure, build, release}`。每次实际改变状态的合法操作，在即将使整个变化持久生效前调用一次 `phase="before_commit"`，在已持久生效但尚未返回前调用一次 `phase="after_commit"`。初次空库初始化、读取、无写入的 build 命中、release 重放不触发。回调异常原样传播（不包装成 `ValueError`）。`before_commit` 处抛出或子进程退出，重开后完全看不到本次变化；`after_commit` 处抛出或子进程退出，重开后完整看到变化。输入错误抛 `ValueError` 且无可见副作用；真实存储故障保留原异常。

## 示例

`examples/archive_sample.py` 在临时目录中演示：建立血缘 A（参考像素 0/10/20/30 → 400/421/444/469，容差 0），登记曝光 EA（`pixel_edges=["0","10","20"]`、`counts=["21","23"]`、`wavelength_edges=["400","410","421","430","444"]`），build 得到精确计数 `["10","11","9","14"]`、`dropped_counts="0"`；随后更正、观察 stale、重建、交付与重放。运行：

    PYTHONPATH=src python3 examples/archive_sample.py

## 范围边界

不处理多谱段、波长反转、峰检测、谱线识别、异常点剔除、并发压测、外部数据库、统计误差估计或大规模性能指标。
