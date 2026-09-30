# PETL：参数化节点增强桁架晶格

PETL 是以下论文的配套研究代码：

> Haochen Li, Bingteng Sun, Lingxin Cao, Xin Yan, and Lin Lu.  
> “A parametric model of joint-enhanced truss lattices with optimization for superior mechanical properties.”  
> *Materials & Design*, 254 (2025), 113969。  
> [论文与 DOI](https://doi.org/10.1016/j.matdes.2025.113969)

[English README](README.md)

PETL 用于生成可调节杆件半径和节点平滑度的桁架晶格，支持隐式几何建模、基于体素的周期均匀化，以及使用移动渐近线法（MMA）的水平集参数优化。模型包含四面体对称（tet）和立方体对称（cubic）。

## 编译

需要 CMake 3.24 或更新版本，以及支持 C++17 的编译器。Windows 请安装 Visual Studio 2019 或更新版本的“使用 C++ 的桌面开发”工作负载和 Windows SDK；Linux 需要 C++ 编译器与 CMake；macOS 需要 Xcode 命令行工具和 CMake。OpenMP 可选。运行集成测试需要 Python 3。

### Windows（PowerShell）

~~~powershell
cmake -S . -B build -A x64
cmake --build build --config Release --parallel
ctest --test-dir build -C Release --output-on-failure
~~~

### Linux / macOS

~~~sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
~~~

CMake 配置会解析 Eigen、libigl 和 MMA 依赖，具体来源与版本见 [cmake/Dependencies.cmake](cmake/Dependencies.cmake)。添加 `-DPETL_ENABLE_OPENMP=OFF` 可关闭 OpenMP，添加 `-DBUILD_TESTING=OFF` 可关闭测试。

## 快速运行

从仓库根目录运行。以下命令使用 Linux / macOS 可执行文件。Windows 请将 `./build/bin/petl` 替换为 `./build/bin/Release/petl.exe`。

生成晶格：

~~~sh
./build/bin/petl generate --skeleton examples/skeletons/bcc.csv --resolution 16 --mesh-resolution 16 --radius 4 --output output/bcc
~~~

计算周期均匀化刚度：

~~~sh
./build/bin/petl homogenize --skeleton examples/skeletons/bcc.csv --resolution 16 --mesh-resolution 16 --radius 4 --output output/homogenized
~~~

运行两次优化迭代：

~~~sh
./build/bin/petl optimize --skeleton examples/skeletons/bcc.csv --family cubic --variant solid --resolution 16 --mesh-resolution 16 --volume-fraction 0.17 --max-iterations 2 --output output/optimization
~~~

运行 `./build/bin/petl --help` 查看全部命令和参数。编译后，Windows 可执行文件位于 `build/bin/Release/petl.exe`；Linux / macOS 可执行文件位于 `build/bin/petl`。

## 输入与参数

骨架文件须包含**完整晶胞桁架**。每条边由连续两行表示，每行包含三个归一化至 [-1, 1] 的坐标。文件不含表头，可使用逗号或空白分隔。坐标必须为有限数值，边长不能为零，输入骨架需符合所选对称类型。examples/skeletons/ 中的文件是用于演示的合成示例。

可选配置文件用于指定 `input_path`、`output_path`、`resolution` 和 `resolution_mesh`。使用 `--config config.txt` 显式加载；配置中的相对路径以配置文件所在目录为基准，命令行参数优先。格式示例见 [config.txt](config.txt)。

| 参数 | 默认值 | 说明 |
| --- | --- | --- |
| `--family` | `tet` | `tet` 四面体对称或 `cubic` 立方体对称 |
| `--variant` | `modified` | `modified` 去除节点鼓包；`solid` 使用原始 KS 构造 |
| `--resolution` | 32 | 每个方向的体素分辨率（4–512） |
| `--mesh-resolution` | 32 | 提取表面网格时的分辨率（4–512） |
| `--radius` | 2.0 | 生成和均匀化使用的统一杆件半径 |
| `--smooth` | 0.5 | 初始节点平滑度 |
| `--volume-fraction` | 生成时不设置；优化时为 0.15 | 目标体积分数；生成或均匀化时据此拟合统一半径 |
| `--max-iterations` | 120 | MMA 最大更新次数，仅用于优化 |

优化器根据 `--family` 和 `--variant` 选择四种预设：

| 对称类型 | 模式 | 优化目标 | 附加约束 |
| --- | --- | --- | --- |
| `tet` | `solid` | 最大化 `C11 + C22 + C33` | 通用各向异性指标 `Au <= 1e-5` |
| `tet` | `modified` | 最大化 `C11 + C22 + C33` | 通用各向异性指标 `Au <= 1e-5` |
| `cubic` | `solid` | 最大化 `C11 + C22 + C33` | 无 |
| `cubic` | `modified` | 最小化 `Au` | 无 |

所有预设均将体积分数限制在 [0.98 * 目标值, 目标值]。半径变量限制在 [0.5, 6]，平滑度限制在 [0.1, 4]。优化可能在达到最大迭代次数前停止。短算例用于检查运行流程，不代表优化已收敛。低分辨率下目标体积分数可能无法达到，半径拟合会失败。

## 输出

`generate` 和 `homogenize` 输出 `mesh.obj`、`voxels.csv` 和 `run.txt`。体素文件包含 `resolution^3` 个二值数据，每行一个值。`homogenize` 还会输出 `stiffness.csv`，即 6×6 等效刚度矩阵；计算使用基体材料参数 E = 1、nu = 0.35，Voigt 顺序为 xx、yy、zz、xy、yz、xz。`optimize` 会输出初始状态和每次迭代的网格、体素文件、`opt.txt` 优化日志及 `run.txt`。

模型采用的晶胞边长为 50，半径和平滑度按该内部长度尺度计量。OBJ 网格会居中并按各坐标轴归一化到 `[-1, 1]`，因此坐标不是物理长度单位。请为不同运行指定不同输出目录；同名文件可能被覆盖。

## 引用

如果在研究中使用 PETL，请引用论文：

~~~bibtex
@article{li2025petl,
  title   = {A parametric model of joint-enhanced truss lattices with optimization for superior mechanical properties},
  author  = {Li, Haochen and Sun, Bingteng and Cao, Lingxin and Yan, Xin and Lu, Lin},
  journal = {Materials \& Design},
  volume  = {254},
  pages   = {113969},
  year    = {2025},
  doi     = {10.1016/j.matdes.2025.113969}
}
~~~

## 许可证

PETL 使用 [MIT 许可证](LICENSE)。第三方软件声明见 [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md)。
