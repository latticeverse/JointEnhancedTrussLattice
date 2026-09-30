#pragma once
#include <Eigen/Core>
#include <cassert>
#include <igl/signed_distance.h>
#include <fstream>
#include <stdexcept>
#include <string>

inline void meshVoxelize(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::VectorXi& ranges, Eigen::MatrixXd& grid, Eigen::VectorXi& label, const std::string& abaqus_path = "")
{

	const auto minV = V.colwise().minCoeff();
	const auto maxV = V.colwise().maxCoeff();
	//Compute grid and corresponding value
	grid.resize(ranges(0) * ranges(1) * ranges(2), 3);

	const auto lerp = [&](const double index, const int axis)-> double {
		return minV(axis) + index / static_cast<double>(ranges(axis)) * (maxV(axis) - minV(axis));
	};
#pragma omp parallel for
	for (int z = 0; z < ranges(2); z++) {
		const double z_point = lerp(static_cast<double>(z) + 0.5, 2);

		for (int y = 0; y < ranges(1); y++) {
			const double y_point = lerp(static_cast<double>(y) + 0.5, 1);

			for (int x = 0; x < ranges(0); x++) {
				const double x_point = lerp(static_cast<double>(x) + 0.5, 0);

				grid.row(x + ranges(0) * (y + ranges(1) * z)) << x_point, y_point, z_point;
			}
		}
	}

	igl::SignedDistanceType type = igl::SIGNED_DISTANCE_TYPE_WINDING_NUMBER;
	Eigen::MatrixXd C, N;
	Eigen::VectorXd S;
	Eigen::VectorXi I;
	igl::signed_distance(grid, V, F, type, S, I, C, N);
	label.setZero(S.rows());

	for (int j = 0; j < S.rows(); j++) label(j) = S(j) >= 0 ? 0 : 1;

	if (abaqus_path.empty()) return;
	std::ofstream outfile(abaqus_path);
	if (!outfile) throw std::runtime_error("Cannot write Abaqus file: " + abaqus_path);
	grid.resize((ranges(0) + 1) * (ranges(1) + 1) * (ranges(2) + 1), 3);
	for (int z = 0; z < ranges(2)+1; z++) {
		const double z_point = lerp(static_cast<double>(z) + 0.5, 2);

		for (int y = 0; y < ranges(1)+1; y++) {
			const double y_point = lerp(static_cast<double>(y) + 0.5, 1);

			for (int x = 0; x < ranges(0) + 1; x++) {
				const double x_point = lerp(static_cast<double>(x) + 0.5, 0);
				int tempid = x + (ranges(0) + 1) * (y + (ranges(1) + 1) * z);
				grid.row(tempid) << x_point, y_point, z_point;
				outfile << tempid + 1 << ",  " << grid.row(tempid)(0) << ",  " << grid.row(tempid)(1) << ",  " << grid.row(tempid)(2) << std::endl;
			}
		}
	}
	outfile << "*ELEMENT,TYPE=C3D8,ELSET=auto2" << std::endl;
	int index = 1;
	for (int z = 0; z < ranges(2); z++) {
		for (int y = 0; y < ranges(1); y++) {
			for (int x = 0; x < ranges(0); x++) {
				int idnumber = x + ranges(0) * (y + ranges(1) * z);
				if (label(idnumber) == 1)
				{
					int leftupid = x + (ranges(0) + 1) * (y + (ranges(1) + 1) * z);
					leftupid += 1;
					//outfile << index << ",     " << leftupid << ",     " << leftupid + 1 << ",     " << leftupid + (ranges(0) + 1) << ",     " << leftupid + (ranges(0) + 1) + 1 << ",     " << leftupid + (ranges(0) + 1) * (ranges(1) + 1) << ",     " << leftupid + (ranges(0) + 1) * (ranges(1) + 1) + 1 << ",     " << leftupid + (ranges(0) + 1) * (ranges(1) + 1) + (ranges(0) + 1) << ",     " << leftupid + (ranges(0) + 1) * (ranges(1) + 1) + (ranges(0) + 1) + 1 << std::endl;
					outfile << index << ",     " << leftupid + 1 << ",     " << leftupid << ",     " << leftupid + (ranges(0) + 1) * (ranges(1) + 1) << ",     " << leftupid + (ranges(0) + 1) * (ranges(1) + 1) + 1 << ",     " << leftupid + (ranges(0) + 1) + 1 << ",     " << leftupid + (ranges(0) + 1) << ",     " << leftupid + (ranges(0) + 1) * (ranges(1) + 1) + (ranges(0) + 1) << ",     " << leftupid + (ranges(0) + 1) * (ranges(1) + 1) + (ranges(0) + 1) + 1 << std::endl;
					index += 1;
				}
			}
		}
	}
	outfile.close();
}
