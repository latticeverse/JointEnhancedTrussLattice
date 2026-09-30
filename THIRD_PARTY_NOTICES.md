# Third-party notices

PETL's project license is the **MIT License**; see [LICENSE](LICENSE). Upstream dependencies retain their own licenses and copyright notices. Their source is downloaded into the ignored build directory, or supplied separately for an offline build.

| Component | Source / revision | License and attribution |
| --- | --- | --- |
| Eigen, including NNLS | [Eigen 5.0.0](https://gitlab.com/libeigen/eigen/-/tree/5.0.0) | Primarily MPL-2.0; some files carry other compatible notices. See upstream `COPYING.README`, `COPYING.MPL2`, and file headers. |
| Eigen NNLS | `unsupported/Eigen/NNLS` in Eigen 5.0.0 | MPL-2.0. Copyright Essex Edwards (2021) and Hannes Matuschek (2013). |
| libigl | [v2.5.0](https://github.com/libigl/libigl/tree/fdaac01bcc52888994f7afd029dcc045dd408484) | Core under MPL-2.0; the `copyleft` module is designated GPL by upstream. See `LICENSE.MPL2`, `LICENSE.GPL`, and file-specific notices. |
| Marching-cubes lookup tables | `include/igl/copyleft/marching_cubes_tables.h` | IsoEx; Computer Graphics Group, RWTH Aachen (2002). The file specifies GNU Library General Public License version 2. Preserve this notice when redistributing it. |
| MMA solver | [jdumas/mma](https://github.com/jdumas/mma/tree/aa51333f942220ac98e5957accb1b7e60590ec6f) | MIT. Copyright Jérémie Dumas (2018). Source credits Niels Aage (original code, 2013), Jun Wu (OpenMP modifications, 2017), and Jérémie Dumas (modifications, 2017). |

The `MMASolver.h` and `MMASolver.cpp` files in the supplied project matched the pinned MMA revision, apart from line endings. They are now fetched rather than copied into the project. OpenGL, GLFW, GLAD, upstream documentation, and precompiled Windows libraries are not part of this release.

The original project CSV helpers credit **Hao Peng (2019-10-23)**. That attribution is retained in the corresponding project headers; CSV parsing was revised for validation and portable compilation. Do not remove contributor attribution when redistributing the project.

For binary distribution, retain applicable notices and provide corresponding source as required by the relevant licenses, including the dependency revisions used to build it. Moving dependencies out of the source tree does not change those obligations. This file is a component inventory; the authoritative terms are the upstream license texts.

The paper is separately published under CC BY-NC 4.0. That article license does not license the software. The original Panetta/ETH datasets and other research assets are not redistributed here; obtain them from their respective owners under their applicable terms.
