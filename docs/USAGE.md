# 使用说明

## 模型与端点

GNEB寻找静态能量面上的路径，所有图像必须共享同一哈密顿量、格点、原子类型和边界条件。开始前分别自由弛豫两个端点，并检查它们是否是需要研究的局部极小值。本版只核验端点力矩残差，没有Hessian极小值认证功能。

使用正常的VAMPIRE `input`、材料文件和UCF描述模型，然后选择 `sim:program=gneb`。温度必须为0。`sim:integrator=llg-heun`可以保留，但GNEB入口不会调用它推进LLG时间。

支持保守双线性交换（包括张量DMI）、原生各向异性以及静态全局/材料外场。非零温度、检查点自动续接、MPI/GPU、偶极、非保守力矩、约束、高阶交换和微磁粗粒化会被拒绝。

## 链文件

`gneb:chain-file`指定一个文本文件，每个非空、非注释行是图像文件名。路径相对于链文件的目录解析。

```text
image_0.txt
image_1.txt
image_2.txt
...
image_8.txt
```

`path`模式要求至少三个图像，包含固定的首尾端点。不要给相邻图像放入完全相同的构型。相邻同一自旋恰好反向时，测地线方向不唯一，应显式补充中间图像，而不是让程序任意选择。

## 图像格式与坐标

```text
GNEB_SPINS_V1
2
0.5 0.5 0.5 0 0 -1
1.5 0.5 0.5 0 0  1
```

首行是格式标识，第二行是自旋数，后面每行为 `x_A y_A z_A sx sy sz`。没有列标题。自旋数应等于VAMPIRE实际生成的磁性格点数；每个格点必须出现一次。

导入器以1e-6 Å分辨率的坐标键匹配内部原子，故行顺序可以改变，但不能只给自旋而省略坐标。使用实际计算坐标；显示、剪切或绘图后的坐标不能直接代替计算坐标。需要更高精度或不规则格点时，应先检查该匹配精度是否能唯一识别所有点。

每个自旋应是单位矢量。文本导出造成的范数偏离≤1e-5时允许单位化并记录最大调整；更大偏离或非有限值拒绝。这个导入不恢复旧RNG状态，不能称为精确轨迹续接。

可选转换工具：

```bash
python3 -m pip install numpy
python3 tools/convert_snapshot.py --coords atoms-coords.data --spins spins-00000050.data --out image.txt
```

此工具针对标准VAMPIRE文本格式：坐标文件跳过首行，取第3–5列；自旋文件跳过首行，取三列自旋。使用前核对实际文件格式与完整性。转换工具只整理数据，不改变哈密顿量，不运行求解器。

## 参数表

所有GNEB参数均使用下面明确的单位，不接受另加 `!unit` 后缀。

|参数|默认值|含义|
|---|---:|---|
|`gneb:mode`|`path`|`path`优化路径；`relax`自由弛豫单个端点；`evaluate`只评估图像|
|`gneb:chain-file`|必须提供|链文件|
|`gneb:threads`|1|OpenMP 原生核线程数；需专门编译，见 ENDPOINTS.md|
|`gneb:checkpoint-every`|0|relax 模式保存构型的迭代间隔；0关闭|
|`gneb:output-prefix`|`gneb`|仅字母、数字、下划线、连字符|
|`gneb:spring-constant`|1|固定弹簧常数，meV/rad²|
|`gneb:force-tolerance`|1e-5|内部像最大GNEB力，meV/rad|
|`gneb:endpoint-tolerance`|1e-5|固定端点最大物理切向力矩，meV/rad|
|`gneb:max-iterations`|20000|有限优化迭代预算|
|`gneb:climbing`|`true`|启用最高能内部像CI|
|`gneb:climb-after`|100|正常预热后启用CI；普通带提前收敛时可提前启用|
|`gneb:report-every`|50|历史记录间隔；最终状态总会记录|
|`gneb:initial-step`|0.05|FIRE初始伪时间步|
|`gneb:maximum-step`|0.25|FIRE最大伪时间步，不小于初始值|
|`gneb:maximum-rotation`|0.1|每次更新单自旋的最大旋转角，rad，范围(0,0.5]|

这些默认值是软件默认，不是所有模型的最佳参数。公开16×16示例使用初始步0.005、最大步0.05，具体见其 `input`。不能把FIRE迭代数或伪时间转换为ns。

## 能量和力约定

程序从VAMPIRE读取的自旋场以Tesla表示，再按各格点磁矩转换为 `−∂E/∂S`，输出meV单位。双向交换表的总能乘1/2，单点项不乘1/2。角度以rad表示，路径长度为

\[
\ell(A,B)=\sqrt{\sum_i\operatorname{atan2}(|\mathbf{s}_{Ai}\times\mathbf{s}_{Bi}|,
\mathbf{s}_{Ai}\cdot\mathbf{s}_{Bi})^2}.
\]

切线、投影和弹簧作用在整幅图像的乘积自旋球面上；它不是逐原子位置的弹簧。普通像的GNEB力由垂直路径的物理切向力和沿路径的弹簧力组成；CI反转沿路径的物理力并去除弹簧。

启动时对确定性选择的格点做切向有限差分，核对原生能量和梯度。这是实现一致性检查，不能取代完整模型与边界条件验证。

## 输出

|文件|含义|
|---|---|
|`*_history.csv`|迭代数、最大GNEB力、端点残差、CI索引|
|`*_energies.csv`|各像测地线反应坐标、总能、相对首像能量、物理切向力矩|
|`*_image_N.txt`|优化后图像，仍为坐标标记格式|
|`*_status.json`|实际数值终态与限制|
|`*_force_N.tsv`|只在静态模式输出逐自旋负Cartesian能量梯度|

`physical_torque_meV`和 `max_GNEB_force_meV`不是同一量：在路径内部，物理力可以沿路径不为零；GNEB力收敛要求其垂直路径分量与图像间距得到控制。

`sampled_forward_barrier_meV`取路径上最高采样能量减首像能量；反向势垒减末像能量。未收敛时这只是末路径的原始值。`evaluate`模式两者均为 `null`。本版不通过能量插值隐含修正势垒。

## 静态评估

```text
gneb:mode=evaluate
gneb:chain-file=one-or-more-images.txt
```

可以只给一个图像。该模式不推进自旋、不优化路径，适合与其他程序交叉核对能量和逐自旋梯度。`numerically_converged`为false，避免把静态评估当成路径收敛。

## 失败与验收

- 退出1：缺文件、模型不支持、非单位/非有限自旋、格点不完整、重复像、反向相邻自旋、端点未驻定、输出已存在等。
- 退出2：耗尽迭代预算，保存实际末路径但 `numerically_converged=false`。
- 路径退出0后仍须检查力、端点、范数、图像数收敛和路径机制。

寻找最低势垒通常还需要不同路径初猜和鞍点负模核验。该扩展不自动生成热寿命、尝试频率或最低路径证明。
