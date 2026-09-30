/* Read an Eigen matrix from CSV.
 * Original helper: Hao Peng, 2019-10-23.
 * Validated parsing and submatrix indexing added during release preparation.
 */
#pragma once
#include <Eigen/Core>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace elib {
// Ranges are zero-based and half-open; -1 selects the corresponding boundary.
template<typename Derived>
bool readCSV(const std::string& filename, int start_row, int start_col,
             int end_row, int end_col, Eigen::PlainObjectBase<Derived>& matrix)
{
    std::ifstream input(filename);
    if (!input) return false;
    std::vector<std::vector<typename Derived::Scalar>> rows;
    std::string line;
    while (std::getline(input, line)) {
        if (line.find_first_not_of(" \t\r") == std::string::npos) continue;
        std::vector<typename Derived::Scalar> row;
        std::istringstream record(line);
        std::string field;
        const bool comma_separated = line.find(',') != std::string::npos;
        while (comma_separated ? static_cast<bool>(std::getline(record, field, ','))
                               : static_cast<bool>(record >> field)) {
            std::istringstream value_stream(field);
            double value;
            if (!(value_stream >> value) || !std::isfinite(value)) return false;
            value_stream >> std::ws;
            if (!value_stream.eof()) return false;
            row.push_back(static_cast<typename Derived::Scalar>(value));
        }
        if (line.back() == ',' || row.empty()) return false;
        if (!rows.empty() && row.size() != rows.front().size()) return false;
        rows.push_back(std::move(row));
    }
    if (rows.empty() || input.bad()) return false;
    const int nrows = static_cast<int>(rows.size());
    const int ncols = static_cast<int>(rows.front().size());
    if (start_row < 0) start_row = 0;
    if (start_col < 0) start_col = 0;
    if (end_row < 0) end_row = nrows;
    if (end_col < 0) end_col = ncols;
    if (start_row >= end_row || start_col >= end_col || end_row > nrows || end_col > ncols)
        return false;
    matrix.resize(end_row - start_row, end_col - start_col);
    for (int r = start_row; r < end_row; ++r)
        for (int c = start_col; c < end_col; ++c)
            matrix(r - start_row, c - start_col) = rows[r][c];
    return true;
}

template<typename Derived>
bool readCSV(const std::string& filename, Eigen::PlainObjectBase<Derived>& matrix)
{
    return readCSV(filename, -1, -1, -1, -1, matrix);
}
} // namespace elib
