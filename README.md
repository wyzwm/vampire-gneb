# VAMPIRE GNEB

为 **VAMPIRE 7.0.0** 提供 `sim:program=gneb` 扩展入口，计算单位自旋系统的最小能量路径和候选跃迁势垒。

本项目参考 [Spirit 的 GNEB 实现](https://github.com/spirit-code/spirit/blob/e82250d3b14411c2c2fa292d143f13e3e111ad8c/core/src/engine/Method_GNEB.cpp)，以独立 C++11 内核实现球面投影、测地线弹簧、改善切线、爬山像和 FIRE 优化，并直接调用 VAMPIRE 原生的能量与有效场代码。**当前版本为单进程 CPU 非官方研究扩展；0.2 支持自由端点弛豫和可选 OpenMP 原子循环并行。**

|入口|内容|
|---|---|
|[安装与编译](docs/INSTALL.md)|Linux / Windows、固定上游版本、构建要求|
|[如何接入扩展入口](docs/INTEGRATION.md)|新增文件、解析路由、程序ID、调度入口、Makefile|
|[完整使用说明](docs/USAGE.md)|参数、图像格式、单位、输出和失败判据|
|[端点弛豫与 CPU 并行](docs/ENDPOINTS.md)|自由端点、OpenMP 支持范围与串并行核验|
|[公开验证算例](docs/VALIDATION.md)|Spirit 官方16×16斯格明子湮灭测试、出处与验收|
|[来源与许可证](NOTICES.md)|VAMPIRE GPL/BSD、Spirit MIT、方法引用|

## 快速开始

依赖：Git、Python 3.8+、支持C++11的GCC和GNU Make。普通构建与公开算例不需要NumPy，也不需要安装Spirit。

```bash
git clone https://github.com/wyzwm/vampire-gneb.git
cd vampire-gneb
python3 tools/fetch_vampire.py
python3 tools/install_extension.py --source build/vampire
python3 tools/build.py --source build/vampire --jobs 2
python3 tools/run_example.py --exe build/bin/vampire-gneb --output runs/first-example
```

Windows将 `python3` 换为 `python`，可执行文件为 `build/bin/vampire-gneb.exe`。编译器配置见[安装说明](docs/INSTALL.md)。仓库如果为私有，需要GitHub账号具有读取权限。

`fetch_vampire.py` 固定到官方标签 `7.0.0`，并核对提交：

```text
525bc27ee44c525aee229570f30f3d4c61d54f66
```

安装器只接受该版本的干净源码，先核对基线文件、补丁和覆盖文件的SHA，再执行接入。不会把新版补丁盲目打到已有修改的源码上。

## 一个完整、公开的例子

仓库自带 [Spirit官方求解器测试](https://github.com/spirit-code/spirit/blob/e82250d3b14411c2c2fa292d143f13e3e111ad8c/core/test/test_solvers.cpp)中的斯格明子湮灭算例。参数来自其 [`solvers.cfg`](https://github.com/spirit-code/spirit/blob/e82250d3b14411c2c2fa292d143f13e3e111ad8c/core/test/input/solvers.cfg)。

|参数|值|
|---|---:|
|格点|16×16×1，简单方格；x、y周期边界|
|交换|最近邻 J=10 meV|
|DMI|最近邻 D=6 meV，Bloch手性+1|
|单自旋磁矩|2 μB|
|外场|Bz=25 T|
|单轴各向异性|K=0|
|路径|9像，包括两个固定端点|
|初态／终态|已弛豫斯格明子／均匀+z态|

Spirit官方测试中的参考总能为：

```text
初态极小值       -5849.69140625 meV
湮灭鞍点         -5811.5244140625 meV
正向势垒            38.1669921875 meV
```

本扩展的复现势垒约 **38.16697 meV**。验收脚本同时核对极小值、鞍点、势垒记账和最大GNEB力，不只检查退出码。官方参考值为float32常量，具体容差及两软件Bohr磁子常数差异见[验证说明](docs/VALIDATION.md)。

预生成图像来自公开模型的独立复现，生成器在 [`benchmarks/spirit_skyrmion.cpp`](benchmarks/spirit_skyrmion.cpp)。原始公开配置、来源提交和SHA也随仓库保留。仓库的验证不依赖未公开材料、外部计算目录或私人轨迹。

## 最小输入

继续使用VAMPIRE的 `input`、`Co.mat`、`file.ucf`，在 `input` 中选择：

```text
sim:program=gneb
sim:temperature=0
gneb:mode=path
gneb:chain-file=chain.txt
gneb:output-prefix=gneb
gneb:spring-constant=1
gneb:force-tolerance=1e-6
gneb:endpoint-tolerance=1e-6
gneb:max-iterations=200000
gneb:climbing=true
gneb:climb-after=2000
gneb:initial-step=0.005
gneb:maximum-step=0.05
gneb:report-every=500
```

`chain.txt`按路径顺序列出图像文件；首尾固定。每个图像包含自旋的实际计算坐标，导入时按坐标匹配，避免内部原子排序改变造成错配。

```text
GNEB_SPINS_V1
1
0.5 0.5 0.5 0 0 1
```

这是格式示意，不是可独立使用的路径。每行依次为 `x y z sx sy sz`，坐标单位Å；自旋单位长度。实际例子有256行/像。

## 输出与成功判据

- `gneb_status.json`：实际模式、迭代数、收敛状态、爬山像和采样势垒。
- `gneb_history.csv`：迭代过程中的最大路径力与端点力矩。
- `gneb_energies.csv`：每像反应坐标、总能、相对能量及物理力矩。
- `gneb_image_*.txt`：最终路径图像。
- `gneb_force_*.tsv`：仅在 `evaluate` 静态模式中输出逐格点负能量梯度。

单位：能量为meV；路径力为meV/rad；反应坐标为乘积自旋球面上的测地线长度。优化器的步长不是LLG物理时间。`evaluate`不会把静态评估写成已收敛路径，其势垒字段为 `null`。

退出码：`0`=成功结束，`1`=输入或模型检查失败，`2`=达到迭代预算但未收敛。对于路径，必须进一步检查 `numerically_converged: true`、力阈值及端点状态。重新运行需新的输出目录或前缀。

## 已实现与适用范围

- 固定端点GNEB、能量加权切线、CI、球面指数映射与速度平行输运。
- VAMPIRE原生双线性交换（包括张量DMI）、各向异性及静态外场。
- 端点残差门禁、坐标/范数核验、相邻重复像和反向自旋拒绝、有限差分一致性检查。
- 单进程CPU。MPI/GPU、偶极、高阶交换、粗粒化、非保守力矩、非零温度和自动检查点续接被明确拒绝。

端点力矩小只证明驻定残差小，不能代替局部极小值/Hessian检查。路径收敛也不能单独证明最低湮灭势垒、鞍点阶数或热寿命。不同初始路径、图像数和体系尺寸仍需针对具体问题检验。

## 测试

```bash
python3 tests/run_tests.py --exe build/bin/vampire-gneb
```

测试包括公开9像与17像基准、解析内核测试、坐标行序置换、非驻定端点、重复/反向图像、非零温度、非单位自旋、迭代上限和输出覆盖保护。GitHub Actions使用同一套公开数据测试。

## 许可证与引用

扩展代码：GPL-2.0-or-later；保留VAMPIRE原有GPL/BSD通知。Spirit参考代码和测试文件为MIT许可，完整许可文本在 [`LICENSES/Spirit-MIT.txt`](LICENSES/Spirit-MIT.txt)。详见 [NOTICES.md](NOTICES.md)。

使用该方法时应引用VAMPIRE、Spirit、GNEB和FIRE的原始方法文献；引用列表及DOI已列在来源说明中。
