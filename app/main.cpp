#include "petl/Truss.h"
#include "petl/homo3D.h"
#include "petl/plsm.h"
#include <cassert>
#include <igl/writeOBJ.h>
#include <filesystem>
#include <iomanip>
#include <map>
#include <set>
#include <stdexcept>

namespace fs = std::filesystem;

namespace {
void usage()
{
    std::cout << R"(PETL - Parametric Joint-Enhanced Truss Lattices
Usage: petl <generate|homogenize|optimize> [options]

  --config FILE          Optional configuration; paths are relative to FILE
  --skeleton FILE        CSV of normalized edge endpoint pairs (required unless
                         --config provides input_path)
  --family tet|cubic     Tetrahedral (1/48) or cubic (1/8) symmetry [tet]
  --variant modified|solid  Enable bulging removal or original KS model [modified]
  --id N                 Dataset id used for input/output names [1]
  --output DIR           Output directory [output], or config output_path
  --resolution N         Voxel resolution, 4..512 [32], or config resolution
  --mesh-resolution N    Mesh resolution, 4..512 [32], or config resolution_mesh
  --smooth S             Joint smoothness, > 0 [0.5]
  --radius R             Initial/uniform radius, > 0 [2.0]
  --volume-fraction V    Fit uniform radius to V, 0 < V < 1;
                         optimization uses target V [0.15]
  --max-iterations N     Optimization update limit [120]
  --help                 Show this help

Examples (run from the repository root):
  petl generate --skeleton examples/skeletons/bcc.csv --output output/bcc
  petl homogenize --skeleton examples/skeletons/bcc.csv --resolution 16 --radius 4
  petl optimize --skeleton examples/skeletons/bcc.csv --family cubic --variant solid --resolution 16 --volume-fraction 0.17 --max-iterations 2

Homogenization uses E = 1 and nu = 0.35. The unit-cell side length is 50.
The included examples are synthetic smoke tests, not the paper's datasets.
)";
}

double number(const std::string& text, const std::string& name)
{
    std::size_t used = 0;
    double value;
    try { value = std::stod(text, &used); }
    catch (...) { throw std::runtime_error("Invalid number for " + name + ": " + text); }
    if (used != text.size() || !std::isfinite(value))
        throw std::runtime_error("Invalid number for " + name + ": " + text);
    return value;
}

int integer(const std::string& text, const std::string& name, int lower, int upper)
{
    const auto value = number(text, name);
    if (value != std::floor(value) || value < lower || value > upper)
        throw std::runtime_error(name + " must be an integer in [" + std::to_string(lower) + ", " + std::to_string(upper) + "]");
    return static_cast<int>(value);
}

void save_geometry(Truss& truss, const fs::path& output)
{
    if (!truss.V.allFinite() || truss.V.rows() == 0 || truss.F.rows() == 0)
        throw std::runtime_error("Empty or non-finite mesh; increase the radius or resolution");
    if (!igl::writeOBJ((output / "mesh.obj").string(), truss.V, truss.F) ||
        !elib::writeCSV((output / "voxels.csv").string(), truss.voxel))
        throw std::runtime_error("Cannot write geometry to " + output.string());
}
} // namespace

int main(int argc, char** argv)
{
    try {
        if (argc == 1 || (argc == 2 && std::string(argv[1]) == "--help")) {
            usage();
            return 0;
        }
        const std::string command = argv[1];
        if (command != "generate" && command != "homogenize" && command != "optimize")
            throw std::runtime_error("Unknown command: " + command);
        const std::set<std::string> allowed = {"--config", "--skeleton", "--family", "--variant", "--id",
            "--output", "--resolution", "--mesh-resolution", "--smooth", "--radius", "--volume-fraction", "--max-iterations"};
        std::map<std::string, std::string> args;
        for (int i = 2; i < argc; ++i) {
            const std::string key = argv[i];
            if (key == "--help") { usage(); return 0; }
            if (!allowed.count(key)) throw std::runtime_error("Unknown option: " + key);
            if (i + 1 == argc || std::string(argv[i + 1]).rfind("--", 0) == 0)
                throw std::runtime_error("Missing value for " + key);
            if (!args.emplace(key, argv[++i]).second) throw std::runtime_error("Duplicate option: " + key);
        }
        auto get = [&](const std::string& key, const std::string& fallback) {
            const auto found = args.find(key);
            return found == args.end() ? fallback : found->second;
        };
        if (command != "optimize" && args.count("--max-iterations"))
            throw std::runtime_error("--max-iterations is only valid for optimize");
        if (command == "optimize" && args.count("--radius"))
            throw std::runtime_error("optimize initializes radius from the target volume fraction; omit --radius");

        Truss truss;
        truss.resolution = truss.resolution_mesh = 32;
        truss.output_path = "output/";
        if (args.count("--config")) truss.loadConfig(args.at("--config"));
        const std::string family = get("--family", "tet");
        const std::string variant = get("--variant", "modified");
        if (family != "tet" && family != "cubic") throw std::runtime_error("--family must be tet or cubic");
        if (variant != "modified" && variant != "solid") throw std::runtime_error("--variant must be modified or solid");
        const bool cubic = family == "cubic", modified = variant == "modified";
        const int id = integer(get("--id", "1"), "--id", 0, 1000000000);
        const int max_iterations = integer(get("--max-iterations", "120"), "--max-iterations", 1, 100000);
        truss.resolution = integer(get("--resolution", std::to_string(truss.resolution)), "--resolution", 4, 512);
        truss.resolution_mesh = integer(get("--mesh-resolution", std::to_string(truss.resolution_mesh)), "--mesh-resolution", 4, 512);
        const double smooth = number(get("--smooth", "0.5"), "--smooth");
        const double radius = number(get("--radius", "2.0"), "--radius");
        const double volume = number(get("--volume-fraction", "0.15"), "--volume-fraction");
        if (smooth <= 0 || radius <= 0) throw std::runtime_error("Smoothness and radius must be positive");
        if (command == "optimize" && (smooth < 0.1 || smooth > 4))
            throw std::runtime_error("Optimization smoothness must be in [0.1, 4]");
        if (volume <= 0 || volume >= 1) throw std::runtime_error("Volume fraction must be in (0, 1)");

        std::string skeleton = get("--skeleton", "");
        if (skeleton.empty()) {
            if (truss.input_path.empty()) throw std::runtime_error("Provide --skeleton or input_path in --config");
            skeleton = (fs::path(truss.input_path) / (std::string(cubic ? "T-E-Skeleton-" : "T-P-Skeleton-") + std::to_string(id) + ".csv")).string();
        }
        if (!fs::is_regular_file(skeleton)) throw std::runtime_error("Skeleton file does not exist: " + skeleton);
        const auto output = fs::absolute(get("--output", truss.output_path)).lexically_normal();
        truss.output_path = output.generic_string() + "/";

        truss.setRInner(0, -1);
        truss.setROuter(radius, -1);
        truss.setSmoothFactor(smooth, -1);
        if (cubic) truss.getSkeletonFromETH(skeleton);
        else truss.getSkeletonFromPanetta(skeleton);
        if (truss.edges.empty() || truss.nodes.empty()) throw std::runtime_error("Skeleton produced no usable topology");
        fs::create_directories(output);

        if (command == "optimize") {
            if (cubic && modified) plsmtoSolidModifiedETH(truss, volume, id, smooth, skeleton, max_iterations);
            else if (cubic) plsmtoSolidETH(truss, volume, id, smooth, skeleton, max_iterations);
            else if (modified) plsmtoSolidModified(truss, volume, id, smooth, skeleton, max_iterations);
            else plsmtoSolid(truss, volume, id, smooth, skeleton, max_iterations);
        } else {
            if (args.count("--volume-fraction")) {
                bool fitted;
                if (cubic && modified) fitted = truss.setSolidVolumeFractionModifiedETH(volume);
                else if (cubic) fitted = truss.setSolidVolumeFractionETH(volume);
                else if (modified) fitted = truss.setSolidVolumeFractionModified(volume);
                else fitted = truss.setSolidVolumeFraction(volume);
                if (!fitted) throw std::runtime_error("Unable to fit the requested volume fraction");
            }
            if (cubic && modified) { truss.genSolidTrussVoxelModifiedETH(); truss.genSolidTrussMeshModifiedETH(); }
            else if (cubic) { truss.genSolidTrussVoxelETH(); truss.genSolidTrussMeshETH(); }
            else if (modified) { truss.genSolidTrussVoxelModified(); truss.genSolidTrussMeshModified(); }
            else { truss.genSolidTrussVoxel(); truss.genSolidTrussMesh(); }
            save_geometry(truss, output);
            if (command == "homogenize") {
                if (truss.voxel.sum() == 0) throw std::runtime_error("No occupied voxels; increase resolution or radius");
                truss.genTrussGrid();
                homo3D::hex_element element;
                homo3D::setMaterial(1.0, 0.35, element);
                homo3D::init(truss.grid, truss.voxel, Eigen::Vector3i::Constant(truss.resolution + 1), element);
                homo3D::solve(element);
                if (!element.macro_constructive_tensor.allFinite()) throw std::runtime_error("Non-finite homogenized tensor");
                if (!elib::writeCSV((output / "stiffness.csv").string(), element.macro_constructive_tensor))
                    throw std::runtime_error("Cannot write stiffness.csv");
                std::cout << "Homogenized stiffness (E=1, nu=0.35):\n" << element.macro_constructive_tensor << '\n';
            }
        }
        std::ofstream report(output / "run.txt");
        report << std::setprecision(17)
               << "command " << command << "\nskeleton " << fs::absolute(skeleton).generic_string()
               << "\nfamily " << family << "\nvariant " << variant << "\nresolution " << truss.resolution
               << "\nmesh_resolution " << truss.resolution_mesh << "\ninitial_radius " << truss.r_outer
               << "\ninitial_smoothness " << smooth << "\nvolume_fraction " << truss.getVolumeFraction()
               << "\ncell_length " << truss.cell_length << '\n';
        report.close();
        if (!report) throw std::runtime_error("Cannot write run.txt");
        std::cout << "Volume fraction: " << truss.getVolumeFraction() << "\nOutput: " << output.string() << '\n';
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "PETL: " << error.what() << '\n';
        return 1;
    }
}
