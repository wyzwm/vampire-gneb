# 安装与编译

## 1. 依赖

- Git；Python 3.8或以上。
- GNU Make；支持C++11的GCC。
- Linux可使用发行版自带的GCC/Make。
- Windows可使用MinGW-w64或便携[w64devkit](https://github.com/skeeto/w64devkit)。把其 `bin` 加入当前终端PATH，或给构建脚本传入编译器和Make的路径。

基本构建、示例和验收只用Python标准库。`tools/convert_snapshot.py`是可选的文本导入辅助工具，需要NumPy。Spirit是算法和公开基准来源，不是运行依赖。

## 2. 获取固定版本

在本仓库根目录执行：

```bash
python3 tools/fetch_vampire.py --destination build/vampire
```

该脚本下载官方 `7.0.0` 标签，强制核对提交 `525bc27ee44c525aee229570f30f3d4c61d54f66`。已有目标目录会被拒绝，不会重置已有工作。

也可以手动获取：

```bash
git -c core.autocrlf=false clone --depth 1 --branch 7.0.0 https://github.com/richard-evans/vampire.git build/vampire
git -C build/vampire rev-parse HEAD
```

## 3. 安装扩展

```bash
python3 tools/install_extension.py --source build/vampire
```

安装顺序为：核对Git提交 → 拒绝已有受跟踪修改 → 核对基线/补丁/覆盖文件SHA → `git apply --check` → 应用接入补丁 → 复制扩展源码和许可证 → 写入本地安装回执。

同一版本重复执行会核对安装文件后直接返回；不同版本或安装后源码被修改会停止。若要开发代码，使用独立的源码副本，重新形成构建和验证记录，不要把旧测试回执当作新版本通过。

## 4. 编译单进程版本

Linux：

```bash
python3 tools/build.py --source build/vampire --jobs 2
```

Windows PowerShell：

```powershell
python tools/build.py --source build/vampire --jobs 2 --compiler C:/tools/w64devkit/bin/g++.exe --make C:/tools/w64devkit/bin/make.exe
```

`C:/tools/...`是示例工具路径，应替换为自己的安装位置。

产物为 `build/bin/vampire-gneb` 或 `build/bin/vampire-gneb.exe`。日志与构建SHA在 `build/`，默认不会提交到Git。脚本只构建串行版本，未使用 `-ffast-math`，不会安装系统服务或覆盖其他程序。

如果自行使用上游Makefile，Linux也可在已安装扩展的源码目录执行：

```bash
make serial -j2 EXECUTABLE=vampire-gneb
```

## 5. 运行示例和验收

```bash
python3 tools/run_example.py --exe build/bin/vampire-gneb --output runs/demo
python3 tests/run_tests.py --exe build/bin/vampire-gneb
```

Windows在上述可执行文件名后加 `.exe`，并使用 `python`。验收需要编译核心测试，可用 `--cxx` 指定GCC路径。

不要把 `make parallel` 或 `mpirun` 当作本版的多核接口：GNEB模块在MPI构建中会拒绝运行。普通上游功能和未来扩展的并行化是不同事项。

## 常见问题

|现象|处理|
|---|---|
|上游版本不匹配|重新取指定7.0.0提交；不跳过校验|
|存在源代码修改|另建干净源码目录，保留现有修改|
|找不到g++/make|配置当前终端PATH或传入完整工具路径|
|输出已存在|使用新的运行目录或输出前缀|
|Linux旧GCC头文件错误|使用仓库附带的兼容补丁和C++11选项|
|尚未得到程序但有日志|读取 `build/build-result.json` 的退出码，不以文件夹存在判断构建成功|

编译成功后还应运行公开验收。不同编译器、硬件或后续源码修改可能改变数值结果，必须记录对应版本。
