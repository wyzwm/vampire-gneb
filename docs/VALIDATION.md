# 公开基准与验收

本仓库的全部公开验收只使用可在线取得的Spirit测试模型和明确写出的解析单元测试，不依赖私人材料参数或轨迹。

## 主基准：Spirit官方16×16斯格明子湮灭

来源仓库：<https://github.com/spirit-code/spirit>，固定提交：

```text
e82250d3b14411c2c2fa292d143f13e3e111ad8c
```

对应原件：

- [`core/test/test_solvers.cpp`](https://github.com/spirit-code/spirit/blob/e82250d3b14411c2c2fa292d143f13e3e111ad8c/core/test/test_solvers.cpp)：初态极小值、鞍点能量和磁化预期。
- [`core/test/input/solvers.cfg`](https://github.com/spirit-code/spirit/blob/e82250d3b14411c2c2fa292d143f13e3e111ad8c/core/test/input/solvers.cfg)：几何、J、D、磁矩、外场及边界。
- [`core/src/utility/Configurations.cpp`](https://github.com/spirit-code/spirit/blob/e82250d3b14411c2c2fa292d143f13e3e111ad8c/core/src/utility/Configurations.cpp)：半径5、相位−90°的初始斯格明子构造。
- [`core/include/utility/Constants.hpp`](https://github.com/spirit-code/spirit/blob/e82250d3b14411c2c2fa292d143f13e3e111ad8c/core/include/utility/Constants.hpp)：单位换算常数。

仓库 `references/spirit/` 保留未修改的测试源码和配置副本及MIT许可证；SHA见 `upstream.json`。

模型为256个单位自旋的简单方格，x/y周期：最近邻J=10 meV，D=6 meV，Bloch手性+1，μ=2 μB，Bz=25 T，K=0，无偶极项。采用唯一正向键写法：

\[
E=-\sum_{\langle ij\rangle}J\mathbf{s}_i\cdot\mathbf{s}_j
  -\sum_{\langle ij\rangle}\mathbf D_{ij}\cdot(\mathbf{s}_i\times\mathbf{s}_j)
  -\mu B_z\sum_i s_{iz}.
\]

沿+x的D矢量为(6,0,0) meV，沿+y为(0,6,0) meV；反向键取反。VAMPIRE使用四个有向键，能量记账包含1/2。

## 图像如何得到

`benchmarks/spirit_skyrmion.cpp`用公开配方生成半径5初始纹理，在同一公开哈密顿量下以球面FIRE自由弛豫端点，再向均匀+z态做球面插值。给定的九个 `image_*.txt` 是这次独立复现生成的图像，**不是冒称来自Spirit仓库的原始轨迹文件**。

重新生成：

```bash
g++ -std=c++11 -O2 benchmarks/spirit_skyrmion.cpp -o build/prepare-public-benchmark
python3 examples/spirit-skyrmion-16x16/generate.py --generator build/prepare-public-benchmark --out build/regenerated-example
```

Windows给输出可执行文件加 `.exe`。生成器输出目录必须不存在。

还可在新的工作目录直接运行该C++程序，额外验证独立公开模型求解器的整条GNEB路径；默认固定9像。这个辅助求解器用于复现模型，主要验收仍运行完整VAMPIRE入口。

## 预期结果与容差

Spirit测试源码写出的总能常量：

|量|参考值|
|---|---:|
|斯格明子极小值能量|−5849.69140625 meV|
|湮灭鞍点能量|−5811.5244140625 meV|
|两者差，即正向势垒|38.1669921875 meV|
|鞍点平均z磁矩|约1.93314 μB/自旋|

这些总能常量由上游测试声明为 `float`。在约5800 meV的绝对能量上，float32的间距约4.88e−4 meV；本仓库对参考极小值和鞍点使用3e−4 meV绝对容差，而不是把源码中的小数位当成双精度参考精度。

另一个已说明差异是Bohr磁子：Spirit该提交使用0.057883817555 meV/T；VAMPIRE使用原生磁矩转换常数再转meV。输入仍保持μ=2 μB、B=25 T，不悄悄调场以拟合参考数值。实际总能和残差完整记录在验收JSON中。

本扩展得到正向势垒约38.16697 meV。9像与17像都需达到最大GNEB力<1e−6 meV/rad，且两者势垒差<1e−5 meV。CI后的鞍点能量、端点能量及记账同时核对。

## 执行测试

```bash
python3 tests/run_tests.py --exe build/bin/vampire-gneb
```

验收项目：

1. 公开Spirit参数下的9像、17像原生VAMPIRE路径。
2. 自旋球面距离、切平面投影、指数更新和平行输运。
3. 额外构造的解析单/双自旋势垒单元测试，公式直接写在 `tests/test_core.cpp`。它们是补充测试，不冒称为Spirit原始算例。
4. 未驻定端点、重复像、反向相邻自旋、错误坐标、非单位自旋和非零温度的拒绝。
5. 迭代预算耗尽的退出2与 `numerically_converged=false`。
6. 输出覆盖拒绝，及坐标行序置换不改变静态能量。

默认验收回执写入被Git忽略的 `build/acceptance.json`。仓库中的 `validation/public-benchmark-results.json` 是经过检查的公开基准数值记录，只包含通用软件/公开模型信息。

## 结果能支持的范围

本测试支持实现能够重现这个公开离散自旋基准和解析测试的数值结果。它不证明任意模型的最低势垒或材料稳定性。不同初猜可能收敛到不同路径；对新的研究问题仍应核对端点局部极小值、鞍点阶数、图像数、边界及尺寸。
