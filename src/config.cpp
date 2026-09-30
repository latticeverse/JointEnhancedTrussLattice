#include "petl/Truss.h"
#include <filesystem>
#include <iomanip>
#include <stdexcept>

bool Truss::loadConfig(std::string filename)
{
    std::ifstream input(filename);
    if (!input) throw std::runtime_error("Cannot open config: " + filename);
    const auto base = std::filesystem::absolute(filename).parent_path();
    std::string line;
    int line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        std::istringstream row(line);
        std::string key;
        if (!(row >> key) || key.front() == '#') continue;
        if (key == "input_path" || key == "output_path") {
            std::string value;
            // Backslashes in a Windows path are path separators, not escapes.
            if (!(row >> std::quoted(value, '"', '\0')) || value.empty())
                throw std::runtime_error("Missing path at config line " + std::to_string(line_number));
            auto path = std::filesystem::path(value);
            if (path.is_relative()) path = base / path;
            (key == "input_path" ? input_path : output_path) = path.lexically_normal().generic_string() + "/";
        } else if (key == "resolution" || key == "resolution_mesh") {
            int value;
            if (!(row >> value) || value < 4 || value > 512)
                throw std::runtime_error("Resolution must be an integer in [4, 512]");
            (key == "resolution" ? resolution : resolution_mesh) = value;
        } else {
            throw std::runtime_error("Unknown config key: " + key);
        }
        std::string extra;
        if ((row >> extra) && extra.front() != '#')
            throw std::runtime_error("Unexpected value at config line " + std::to_string(line_number));
    }
    if (input.bad()) throw std::runtime_error("Cannot read config: " + filename);
    return true;
}
