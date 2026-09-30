# PETL

**Parametric Joint-Enhanced Truss Lattices** is the research code accompanying:

> Haochen Li, Bingteng Sun, Lingxin Cao, Xin Yan, and Lin Lu.  
> “A parametric model of joint-enhanced truss lattices with optimization for superior mechanical properties.”  
> *Materials & Design*, 254 (2025), 113969.  
> [Paper and DOI](https://doi.org/10.1016/j.matdes.2025.113969)

[中文说明 / Chinese README](README.zh-CN.md)

PETL generates truss lattices with controllable strut radii and joint smoothness. It supports implicit geometry generation, voxel-based periodic homogenization, and level-set parameter optimization using the Method of Moving Asymptotes (MMA). The model includes tetrahedral (tet) and cubic (cubic) symmetry, with optional removal of joint bulging regions.

## Build

You need CMake 3.24 or newer and a C++17 compiler. On Windows, install Visual Studio 2019 or newer with the **Desktop development with C++** workload and a Windows SDK. On Linux, install a C++ compiler and CMake. macOS requires Xcode Command Line Tools and CMake. The first CMake configuration downloads the pinned dependencies, so an internet connection is required. OpenMP is optional. Python 3 is required to run the integration tests.

### Windows (PowerShell)

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

The CMake configuration resolves the required Eigen, libigl, and MMA dependencies. Their source locations and versions are recorded in [cmake/Dependencies.cmake](cmake/Dependencies.cmake). To disable OpenMP, configure with `-DPETL_ENABLE_OPENMP=OFF`. To disable tests, use `-DBUILD_TESTING=OFF`.

## Quick start

Run commands from the repository root. The examples below use the Linux/macOS executable. On Windows, replace `./build/bin/petl` with `./build/bin/Release/petl.exe`.

Generate a lattice:

~~~sh
./build/bin/petl generate --skeleton examples/skeletons/bcc.csv --resolution 16 --mesh-resolution 16 --radius 4 --output output/bcc
~~~

Compute the periodic homogenized stiffness tensor:

~~~sh
./build/bin/petl homogenize --skeleton examples/skeletons/bcc.csv --resolution 16 --mesh-resolution 16 --radius 4 --output output/homogenized
~~~

Run a short optimization:

~~~sh
./build/bin/petl optimize --skeleton examples/skeletons/bcc.csv --family cubic --variant solid --resolution 16 --mesh-resolution 16 --volume-fraction 0.17 --max-iterations 2 --output output/optimization
~~~

Run `./build/bin/petl --help` for the full command and option list. The executable is written to `build/bin/Release/petl.exe` on Windows and `build/bin/petl` on Linux or macOS.

## Input and parameters

A skeleton file contains the full unit-cell truss. Each edge is specified by two consecutive rows, with three normalized coordinates in [-1, 1] per row. Files have no header and may use comma or whitespace separators. Coordinates must be finite, edges must have nonzero length, and the skeleton must match the selected symmetry family. The files in examples/skeletons/ are small synthetic examples intended to demonstrate the program.

An optional configuration file sets `input_path`, `output_path`, `resolution`, and `resolution_mesh`. Pass it explicitly with `--config config.txt`; relative paths in the configuration are resolved from the configuration file. Command-line options override configured values. See [config.txt](config.txt) for an example.

| Option | Default | Description |
| --- | --- | --- |
| `--family` | `tet` | `tet` (tetrahedral symmetry) or `cubic` |
| `--variant` | `modified` | `modified` removes joint bulging; `solid` uses the original KS construction |
| `--resolution` | 32 | Voxel resolution per axis (4–512) |
| `--mesh-resolution` | 32 | Resolution used to extract the surface mesh (4–512) |
| `--radius` | 2.0 | Uniform strut radius for geometry generation and homogenization |
| `--smooth` | 0.5 | Initial joint smoothness |
| `--volume-fraction` | unset for generation; 0.15 for optimization | Target volume fraction; for generation and homogenization, fits the uniform radius |
| `--max-iterations` | 120 | Maximum MMA updates; optimization only |

The optimizers provide four presets selected by `--family` and `--variant`:

| Family | Variant | Objective | Additional constraint |
| --- | --- | --- | --- |
| `tet` | `solid` | Maximize `C11 + C22 + C33` | Universal anisotropy index `Au <= 1e-5` |
| `tet` | `modified` | Maximize `C11 + C22 + C33` | Universal anisotropy index `Au <= 1e-5` |
| `cubic` | `solid` | Maximize `C11 + C22 + C33` | None |
| `cubic` | `modified` | Minimize `Au` | None |

All presets constrain the volume fraction to [0.98 * target, target]. Optimization variables are bounded to radii [0.5, 6] and smoothness [0.1, 4]. The optimizer may stop before reaching the requested iteration limit. A short run checks the workflow but does not establish convergence. Volume fitting may fail at coarse resolutions when the requested fraction is not attainable.

## Output

`generate` and `homogenize` write `mesh.obj`, `voxels.csv`, and `run.txt`. The voxel file contains `resolution^3` binary values, one per row. `homogenize` also writes `stiffness.csv`, a 6×6 effective stiffness tensor calculated with base material E = 1 and nu = 0.35 (Voigt order: xx, yy, zz, xy, yz, xz). `optimize` writes the initial and per-iteration meshes and voxel files, an `opt.txt` log, and `run.txt` in the requested output directory.

The unit-cell side length used by the model is 50; radii and smoothness use this internal length scale. OBJ meshes are centered and normalized per axis to `[-1, 1]`, so their coordinates are not physical units. Use a separate output directory for each run, since files with matching names may be overwritten.

## Citation

If you use PETL in your research, please cite the paper:

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

## License

PETL is released under the [MIT License](LICENSE). See [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) for third-party software notices.
