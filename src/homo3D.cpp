#include "petl/homo3D.h"
#include <chrono>
void homo3D::init(const Eigen::MatrixXd& grid, const Eigen::VectorXi& labels, const Eigen::Vector3i& ranges, hex_element& hexa)
{
	Eigen::MatrixXd grids = grid;
	int sx = ranges(0), sy = ranges(1), sz = ranges(2);

	/*Assign the node indices of front bound to back bound at three axis

	[1 4 7]	[10 13 16] [19 22 25]		[1 3 1]	[5 7 5] [1 3 1]		y
	[2 5 8] [11 14 17] [20 23 26]  ->	[2 4 2] [6 8 6] [2 4 2]     |
	[3 6 9] [12 15 18] [21 24 27]		[1 3 1] [5 7 5] [1 3 1]     |__ __x

	*/
	std::vector<std::vector<std::vector<int>>> node_indices(sx, std::vector<std::vector<int>>(sy, std::vector<int>(sz, 0)));

	for (int z = 0; z < sz; ++z)
		for (int y = 0; y < sy; ++y)
			for (int x = 0; x < sx; ++x)
				node_indices[x][y][z] = x + sx * (y + sy * z);
	//periodic bound for X-axis
	for (int z = 0; z < sz; ++z)
		for (int y = 0; y < sy; ++y)
			node_indices[sx - 1][y][z] = node_indices[0][y][z];

	//periodic bound for Y-axis
	for (int z = 0; z < sz; ++z)
		for (int x = 0; x < sx; ++x)
			node_indices[x][sy - 1][z] = node_indices[x][0][z];

	//periodic bound for Z-axis
	for (int y = 0; y < sy; ++y)
		for (int x = 0; x < sx; ++x)
			node_indices[x][y][sz - 1] = node_indices[x][y][0];

	// Vertex order for cubic meshes
	//
	//     4 - - - 7
	//    /|      /|
	//   5 - - - 6 |       z
	//   | |     | |       |
	//   | 0 - - | 3       |_ _ _y
	//   |/      |/       /
	//   1 - - - 2       x
	Eigen::MatrixXi corner(8, 3);
	corner << 0, 0, 0,
		1, 0, 0,
		1, 1, 0,
		0, 1, 0,
		0, 0, 1,
		1, 0, 1,
		1, 1, 1,
		0, 1, 1;

	const auto AperiodicIndex = [&sx, &sy, &sz, &corner](const Eigen::Vector3i& cell)->Eigen::RowVectorXi {
		Eigen::MatrixXi member = cell.replicate(1, 8) + corner.transpose();
		return member.row(0) + sx * (member.row(1) + sy * member.row(2));
	};

	const auto PeriodicIndex = [&node_indices, &corner](const Eigen::Vector3i& cell, const int i)->int {
		return node_indices[cell(0) + corner(i, 0)][cell(1) + corner(i, 1)][cell(2) + corner(i, 2)];
	};

	//extract solid voxel (x,y,z) x(0,range) y(0,range) z(0,range)
	std::vector<Eigen::Vector3i> cell_pos;

	std::set<int> aperiodic_indices, periodic_indices;

	//accoording z-y-x order save element index whose voxel==1
	for (int z = 0; z < sz - 1; ++z) {
		for (int y = 0; y < sy - 1; ++y) {
			for (int x = 0; x < sx - 1; ++x) {
				if (labels(x + (sx - 1) * (y + (sy - 1) * z)) == 1) {
					cell_pos.emplace_back(x, y, z);
					auto member = AperiodicIndex(cell_pos.back());
					for (int i = 0; i < 8; ++i) {
						aperiodic_indices.insert(member(i));
						periodic_indices.insert(PeriodicIndex(cell_pos.back(), i));
					}
				}

			}
		}
	}

	// fill nodes
	Eigen::MatrixXd nodes(aperiodic_indices.size(), 3);

	//grid(range*range*range,3) {x,y,z} of voxel center
	//size: number of voxel range*range*range
	//IC save the index to find a vertex in nodes(coordinates of each vertex) IC save the index relation between aperiodic_indices and all indices
	std::vector<int> IC(grids.rows(), -1), periodic_IC(grids.rows(), -1);

	//aperiodic_indeces vertex index of element which voxel==1 size 8*size(voxel==1)b
	int index = 0;
	for (auto it = aperiodic_indices.begin(); it != aperiodic_indices.end(); ++it) {
		nodes.row(index) = grids.row(*it);
		IC[*it] = index++;
	}

	//periodic_IC saves the index relation between periodic_indices and all indices
	index = 0;
	for (auto it = periodic_indices.begin(); it != periodic_indices.end(); ++it) {
		periodic_IC[*it] = index++;
	}

	//fill cells and periodic_cells
	Eigen::MatrixXi cells(cell_pos.size(), 8);//element of voxel==1
	hexa.periodic_cell_dofs.resize(cell_pos.size(), 24); //element of voxel==1, 8*3

	//std::cout << "cellpos_size: " << cell_pos.size() << std::endl;

	for (int i = 0; i < cell_pos.size(); ++i) {
		auto cell = AperiodicIndex(cell_pos[i]);
		for (int j = 0; j < 8; ++j) {
			cells(i, j) = IC[cell(j)];
			for (int k = 0; k < 3; ++k) {
				hexa.periodic_cell_dofs(i, j * 3 + k) = periodic_IC[PeriodicIndex(cell_pos[i], j)] * 3 + k;
			}
		}
	}

	hexa.nodes_ = nodes; //coordinates
	hexa.cells_ = cells; //indices of each 
}

void homo3D::init_simp(const Eigen::MatrixXd& grid, const Eigen::VectorXi& labels, const Eigen::Vector3i& ranges, hex_element& hexa)
{
	int sx = ranges(0), sy = ranges(1), sz = ranges(2);

	/*Assign the node indices of front bound to back bound at three axis

	[1 4 7]	[10 13 16] [19 22 25]		[1 3 1]	[5 7 5] [1 3 1]		y
	[2 5 8] [11 14 17] [20 23 26]  ->	[2 4 2] [6 8 6] [2 4 2]     |
	[3 6 9] [12 15 18] [21 24 27]		[1 3 1] [5 7 5] [1 3 1]     |__ __x

	*/
	std::vector<std::vector<std::vector<int>>> node_indices(sx, std::vector<std::vector<int>>(sy, std::vector<int>(sz, 0)));

	for (int z = 0; z < sz; ++z)
		for (int y = 0; y < sy; ++y)
			for (int x = 0; x < sx; ++x)
				node_indices[x][y][z] = x + sx * (y + sy * z);
	//periodic bound for X-axis
	for (int z = 0; z < sz; ++z)
		for (int y = 0; y < sy; ++y)
			node_indices[sx - 1][y][z] = node_indices[0][y][z];

	//periodic bound for Y-axis
	for (int z = 0; z < sz; ++z)
		for (int x = 0; x < sx; ++x)
			node_indices[x][sy - 1][z] = node_indices[x][0][z];

	//periodic bound for Z-axis
	for (int y = 0; y < sy; ++y)
		for (int x = 0; x < sx; ++x)
			node_indices[x][y][sz - 1] = node_indices[x][y][0];

	// Vertex order for cubic meshes
	//
	//     4 - - - 7
	//    /|      /|
	//   5 - - - 6 |       z
	//   | |     | |       |
	//   | 0 - - | 3       |_ _ _y
	//   |/      |/       /
	//   1 - - - 2       x
	Eigen::MatrixXi corner(8, 3);
	corner << 0, 0, 0,
		1, 0, 0,
		1, 1, 0,
		0, 1, 0,
		0, 0, 1,
		1, 0, 1,
		1, 1, 1,
		0, 1, 1;

	const auto AperiodicIndex = [&sx, &sy, &sz, &corner](const Eigen::Vector3i& cell)->Eigen::RowVectorXi {
		Eigen::MatrixXi member = cell.replicate(1, 8) + corner.transpose();
		return member.row(0) + sx * (member.row(1) + sy * member.row(2));
	};

	const auto PeriodicIndex = [&node_indices, &corner](const Eigen::Vector3i& cell, const int i)->int {
		return node_indices[cell(0) + corner(i, 0)][cell(1) + corner(i, 1)][cell(2) + corner(i, 2)];
	};

	//extract solid voxel (x,y,z) x(0,range) y(0,range) z(0,range)
	std::vector<Eigen::Vector3i> cell_pos;

	std::set<int> aperiodic_indices, periodic_indices;

	Eigen::VectorXd unitv = Eigen::VectorXd::Ones(labels.rows());
	Eigen::VectorXd inc = -1e-3 * (labels.cast<double>() - unitv);
	hexa.densities = labels.cast<double>() + inc;

	std::vector<Eigen::Triplet<double>> triplets;

	int nz = sz - 1;
	int ny = sy - 1;
	int nx = sx - 1;
	for (int z = 0; z < nz; ++z) {
		for (int y = 0; y < ny; ++y) {
			for (int x = 0; x < nx; ++x) {
				cell_pos.emplace_back(x, y, z);
				auto member = AperiodicIndex(cell_pos.back());
				for (int i = 0; i < 8; ++i) {
					aperiodic_indices.insert(member(i));
					periodic_indices.insert(PeriodicIndex(cell_pos.back(), i));
				}

				int e1 = z * nx * ny + y * nx + x;

				for (int z1 = std::max(z - (ceil(hexa.rmin) - 1), (double)0); z1 <= std::min(z + (ceil(hexa.rmin) - 1), (double)(nz - 1)); z1++)
				{
					for (int y1 = std::max(y - (ceil(hexa.rmin) - 1), (double)0); y1 <= std::min(y + (ceil(hexa.rmin) - 1), (double)(ny - 1)); y1++)
					{
						for (int x1 = std::max(x - (ceil(hexa.rmin) - 1), (double)0); x1 <= std::min(x + (ceil(hexa.rmin) - 1), (double)(nx - 1)); x1++)
						{
							int e2 = z1 * nx * ny + y1 * nx + x1;
							triplets.emplace_back(e1, e2, std::max((double)0, hexa.rmin - std::sqrt(pow(x - x1, 2) + pow(y - y1, 2) + pow(z - z1, 2))));
						}
					}
				}
			}
		}
	}

	hexa.H.resize(labels.rows(), labels.rows());

	hexa.H.setFromTriplets(triplets.begin(), triplets.end());
	hexa.H.makeCompressed();

	hexa.Hs.setZero(labels.rows());
#pragma omp parallel for
	for (int i = 0; i < labels.rows(); i++) {
		hexa.Hs(i) = hexa.H.row(i).sum();
	}
	// fill nodes
	Eigen::MatrixXd nodes(aperiodic_indices.size(), 3);

	//grid(range*range*range,3) {x,y,z} of voxel center
	//size: number of voxel range*range*range
	//IC save the index to find a vertex in nodes(coordinates of each vertex) IC save the index relation between aperiodic_indices and all indices
	std::vector<int> IC(grid.rows(), -1), periodic_IC(grid.rows(), -1);

	//aperiodic_indeces vertex index of element which voxel==1 size 8*size(voxel==1)b
	int index = 0;
	for (auto it = aperiodic_indices.begin(); it != aperiodic_indices.end(); ++it) {
		nodes.row(index) = grid.row(*it);
		IC[*it] = index++;
	}

	//periodic_IC saves the index relation between periodic_indices and all indices
	index = 0;
	for (auto it = periodic_indices.begin(); it != periodic_indices.end(); ++it) {
		periodic_IC[*it] = index++;
	}

	//fill cells and periodic_cells
	Eigen::MatrixXi cells(cell_pos.size(), 8);//element of voxel==1
	hexa.periodic_cell_dofs.resize(cell_pos.size(), 24); //element of voxel==1, 8*3


	for (int i = 0; i < cell_pos.size(); ++i) {
		auto cell = AperiodicIndex(cell_pos[i]);
		for (int j = 0; j < 8; ++j) {
			cells(i, j) = IC[cell(j)];
			for (int k = 0; k < 3; ++k) {
				hexa.periodic_cell_dofs(i, j * 3 + k) = periodic_IC[PeriodicIndex(cell_pos[i], j)] * 3 + k;
			}
		}
	}

	hexa.nodes_ = nodes; //coordinates
	hexa.cells_ = cells; //indices of each 
}

void homo3D::setMaterial(const double E, const double v, hex_element& hexa)
{
	double lambda = v / (1 + v) / (1 - 2 * v) * E, mu = 1. / (2. * (1. + v)) * E;
	hexa.micro_constructive_tensor.resize(6, 6);

	//micro_constructive_tensor <<
	//	3, 1, 1, 0, 0, 0,
	//	1, 3, 1, 0, 0, 0,
	//	1, 1, 3, 0, 0, 0,
	//	0, 0, 0, 1, 0, 0,
	//	0, 0, 0, 0, 1, 0,
	//	0, 0, 0, 0, 0, 1;
	hexa.micro_constructive_tensor <<
		lambda + 2 * mu, lambda, lambda, 0, 0, 0,
		lambda, lambda + 2 * mu, lambda, 0, 0, 0,
		lambda, lambda, lambda + 2 * mu, 0, 0, 0,
		0, 0, 0, mu, 0, 0,
		0, 0, 0, 0, mu, 0,
		0, 0, 0, 0, 0, mu;
}

void homo3D::setMaterial_simp(const double v, hex_element& hexa)
{
	double unit_lambda = v / (1 + v) / (1 - 2 * v); double unit_mu = 1. / (2. * (1. + v));

	hexa.unit_mic_c.resize(6, 6);

	hexa.unit_mic_c <<
		unit_lambda + 2 * unit_mu, unit_lambda, unit_lambda, 0, 0, 0,
		unit_lambda, unit_lambda + 2 * unit_mu, unit_lambda, 0, 0, 0,
		unit_lambda, unit_lambda, unit_lambda + 2 * unit_mu, 0, 0, 0,
		0, 0, 0, unit_mu, 0, 0,
		0, 0, 0, 0, unit_mu, 0,
		0, 0, 0, 0, 0, unit_mu;
}

void homo3D::getCellStiffnessAndForce(const int index, const Eigen::MatrixXd& construtive_tensor, Eigen::MatrixXd& stiffness, Eigen::MatrixXd& force, hex_element& hexa)
{
	stiffness.setZero(24, 24);
	force.setZero(24, 6);

	//Gauss numerical integration with 3 points
	Eigen::Vector3d integral_points(-0.774596669241483, 0, 0.774596669241483);
	Eigen::Vector3d weights(0.555555555555556, 0.888888888888889, 0.555555555555556);

	Eigen::MatrixXd q(3, 8), corner(8, 3), jacobian(3, 3), B(6, 24);
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			for (int k = 0; k < 3; ++k) {

				double x = integral_points(i), y = integral_points(j), z = integral_points(k);

				q << -((y - 1) * (z - 1)) / 8, ((y - 1) * (z - 1)) / 8,
					-((y + 1) * (z - 1)) / 8, ((y + 1) * (z - 1)) / 8,
					((y - 1) * (z + 1)) / 8, -((y - 1) * (z + 1)) / 8,
					((y + 1) * (z + 1)) / 8, -((y + 1) * (z + 1)) / 8,
					-((x - 1) * (z - 1)) / 8, ((x + 1) * (z - 1)) / 8,
					-((x + 1) * (z - 1)) / 8, ((x - 1) * (z - 1)) / 8,
					((x - 1) * (z + 1)) / 8, -((x + 1) * (z + 1)) / 8,
					((x + 1) * (z + 1)) / 8, -((x - 1) * (z + 1)) / 8,
					-((x - 1) * (y - 1)) / 8, ((x + 1) * (y - 1)) / 8,
					-((x + 1) * (y + 1)) / 8, ((x - 1) * (y + 1)) / 8,
					((x - 1) * (y - 1)) / 8, -((x + 1) * (y - 1)) / 8,
					((x + 1) * (y + 1)) / 8, -((x - 1) * (y + 1)) / 8;

				for (int m = 0; m < 8; ++m) {
					corner.row(m) = hexa.nodes_.row(hexa.cells_(index, m));
				}

				jacobian.noalias() = q * corner;

				Eigen::MatrixXd qxyz = jacobian.lu().solve(q);

				for (int m = 0; m < 8; m++) {
					Eigen::MatrixXd Be(6, 3);
					Be << qxyz(0, m), 0, 0,
						0, qxyz(1, m), 0,
						0, 0, qxyz(2, m),
						qxyz(1, m), qxyz(0, m), 0,
						0, qxyz(2, m), qxyz(1, m),
						qxyz(2, m), 0, qxyz(0, m);
					B.block(0, 0 + 3 * m, 6, 3) = Be;
				}

				double weight = jacobian.determinant() * weights(i) * weights(j) * weights(k);

				stiffness += weight * B.transpose() * construtive_tensor * B;
				force += weight * B.transpose() * construtive_tensor;
			}
		}
	}
	/*for (int i = 0; i < force.rows(); i++)
	{
		for (int j = 0; j < force.cols(); j++)
			std::cout << force(i,j) << " ";
		std::cout << std::endl;
	}*/
}

void homo3D::get_sigma_f(const int index, Eigen::VectorXd& sigma, Eigen::MatrixXd& f, Eigen::MatrixXd& d, Eigen::VectorXd& m1, hex_element& hexa)
{
	double x = 0.0, y = 0.0, z = 0.0;
	Eigen::MatrixXd q(3, 8), corner(8, 3), jacobian(3, 3), B(6, 24);
	q << -((y - 1) * (z - 1)) / 8, ((y - 1) * (z - 1)) / 8,
		-((y + 1) * (z - 1)) / 8, ((y + 1) * (z - 1)) / 8,
		((y - 1) * (z + 1)) / 8, -((y - 1) * (z + 1)) / 8,
		((y + 1) * (z + 1)) / 8, -((y + 1) * (z + 1)) / 8,
		-((x - 1) * (z - 1)) / 8, ((x + 1) * (z - 1)) / 8,
		-((x + 1) * (z - 1)) / 8, ((x - 1) * (z - 1)) / 8,
		((x - 1) * (z + 1)) / 8, -((x + 1) * (z + 1)) / 8,
		((x + 1) * (z + 1)) / 8, -((x - 1) * (z + 1)) / 8,
		-((x - 1) * (y - 1)) / 8, ((x + 1) * (y - 1)) / 8,
		-((x + 1) * (y + 1)) / 8, ((x - 1) * (y + 1)) / 8,
		((x - 1) * (y - 1)) / 8, -((x + 1) * (y - 1)) / 8,
		((x + 1) * (y + 1)) / 8, -((x - 1) * (y + 1)) / 8;

	for (int m = 0; m < 8; ++m) {
		corner.row(m) = hexa.nodes_.row(hexa.cells_(index, m));
	}

	jacobian.noalias() = q * corner;

	Eigen::MatrixXd qxyz = jacobian.lu().solve(q);

	for (int m = 0; m < 8; m++) {
		Eigen::MatrixXd Be(6, 3);
		Be << qxyz(0, m), 0, 0,
			0, qxyz(1, m), 0,
			0, 0, qxyz(2, m),
			qxyz(1, m), qxyz(0, m), 0,
			0, qxyz(2, m), qxyz(1, m),
			qxyz(2, m), 0, qxyz(0, m);
		B.block(0, 0 + 3 * m, 6, 3) = Be;
	}

	Eigen::MatrixXd f_strain = B * d;

	sigma = hexa.micro_constructive_tensor * f_strain * hexa.macro_compliance_tensor * m1;
	f = hexa.micro_constructive_tensor * f_strain * hexa.macro_compliance_tensor;
}

void homo3D::solve(hex_element& hexa)
{
	clock_t start, end;

	int node_size = hexa.periodic_cell_dofs.maxCoeff() + 1;
	int cell_size = hexa.cells_.rows();

	Eigen::MatrixXd Ke, Fe;
	getCellStiffnessAndForce(0, hexa.micro_constructive_tensor, Ke, Fe, hexa);

	//Assembly
	Eigen::SparseMatrix<double> K(node_size - 3, node_size - 3);
	Eigen::MatrixXd F = Eigen::MatrixXd::Zero(node_size - 3, 6);
	{
		start = clock();

		std::vector<Eigen::Triplet<double>> triplets;

		for (int i = 0; i < cell_size; ++i) {

			for (int r = 0; r < 24; ++r) {
				int index_r = hexa.periodic_cell_dofs(i, r) - 3;
				if (index_r < 0) continue;

				for (int c = 0; c < 24; ++c) {
					int index_c = hexa.periodic_cell_dofs(i, c) - 3;
					if (index_c < 0) continue;
					if (Ke(r, c) != 0)
						triplets.emplace_back(index_r, index_c, Ke(r, c));
				}
				F.row(index_r) += Fe.row(r);
			}
		}

		K.setFromTriplets(triplets.begin(), triplets.end());
		K.makeCompressed();
		//Eigen::MatrixXd K_ = K;
		//std::cout << K.block(0, 0, 50, 50) << std::endl;

		end = clock();

		//std::cout << "Assemble time: " << end - start << " ms" << std::endl;
	}

	//solve KX=f
	Eigen::MatrixXd X;
	{
		start = clock();
		Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Upper> cg;
		//cg.setMaxIterations(200);
		cg.setTolerance(1e-6);
		cg.compute(K);
		X = cg.solve(F);
		end = clock();
		//std::cout << "Solve time: " << (end - start) << " ms" << std::endl;
	}

	start = clock();

	//fix dof= {0,1,2,4,5,11} for cell stiffness, however I do not know why...
	std::vector<int> unfixed_dofs = { 3,6,7,8,9,10,12,13,14,15,16,17,18,19,20,21,22,23 };

	Eigen::MatrixXd Ke_(unfixed_dofs.size(), unfixed_dofs.size()), Fe_(unfixed_dofs.size(), 6);

	for (int r = 0; r < unfixed_dofs.size(); ++r) {
		for (int c = 0; c < unfixed_dofs.size(); ++c) {
			Ke_(r, c) = Ke(unfixed_dofs[r], unfixed_dofs[c]);
		}
		Fe_.row(r) = Fe.row(unfixed_dofs[r]);
	}

	Eigen::MatrixXd temp = Ke_.lu().solve(Fe_);

	Eigen::MatrixXd X0 = Eigen::MatrixXd::Zero(24, 6);
	for (int i = 0; i < unfixed_dofs.size(); ++i)
		X0.row(unfixed_dofs[i]) = temp.row(i);

	//pre-computing
	std::vector<Eigen::MatrixXd> TX0(6), TX(6);
	for (int i = 0; i < 6; i++) {
		TX0[i] = X0.col(i).transpose().replicate(cell_size, 1);

		TX[i].resize(cell_size, 24);

		for (int r = 0; r < cell_size; r++) {
			for (int c = 0; c < 24; c++) {
				TX[i](r, c) = hexa.periodic_cell_dofs(r, c) - 3 < 0 ? 0 : X(hexa.periodic_cell_dofs(r, c) - 3, i);
			}
		}
	}

	hexa.macro_constructive_tensor = Eigen::MatrixXd::Zero(6, 6);


	/*!< \f[C^H_{ij}=1/|\Omega|\sum_e\int_{\Omega_e}(\Chi_e^{0(i)}-\Chi_e^(i))^T k_e (\Chi_e^{0(j)}-\Chi_e^(j))d\Omega^(e)\f] */
	Eigen::MatrixXd left, right;
	for (int i = 0; i < 6; i++) {
		left.noalias() = (TX0[i] - TX[i]) * Ke;

		for (int j = 0; j < 6; j++) {
			right.noalias() = TX0[j] - TX[j];
			Eigen::VectorXd sumrows = (left.array() * right.array()).rowwise().sum();

			hexa.macro_constructive_tensor(i, j) = (left.array() * right.array()).sum();
			if (abs(hexa.macro_constructive_tensor(i, j)) < 1e-5) hexa.macro_constructive_tensor(i, j) = 0;

		}
	}


	hexa.macro_constructive_tensor *= 1 / (hexa.nodes_.colwise().maxCoeff() - hexa.nodes_.colwise().minCoeff()).prod();


	hexa.macro_compliance_tensor = hexa.macro_constructive_tensor.inverse();

	//calculate vms
	//Eigen::MatrixXd d(24, 6);
	//Eigen::MatrixXd f;
	//Eigen::VectorXd sigma;
	////Eigen::VectorXd sigma2(6); sigma2 << 0, 0, 1, 0, 0, 0;
	////hexa.vms_vox.resize(cell_size);

	//for (int i = 0; i < cell_size; i++)
	//{
	//	for (int j = 0; j < 24; j++)
	//	{
	//		for (int r = 0; r < 6; r++) {
	//			d(j, r) = hexa.periodic_cell_dofs(i, j) - 3 < 0 ? 0 : X(hexa.periodic_cell_dofs(i, j) - 3, r);
	//		}
	//	 }

	//	Eigen::MatrixXd td = X0 - d;
	//	//std::cout << d << std::endl;
	//	////std::cout << hexa.micro_constructive_tensor << std::endl;
	//	//system("pause");

	//	//get_sigma_f(0, sigma, f, td, sigma2, hexa);
	//	//std::cout << sigma << std::endl;
	//	//system("pause");

	//	//double vms = std::sqrt(0.5 * (std::pow((sigma(0, 0) - sigma(1, 0)), 2) + std::pow((sigma(1, 0) - sigma(2, 0)), 2) 
	//	//	+ std::pow((sigma(2, 0) - sigma(0, 0)), 2)) + 3 * (std::pow(sigma(3, 0), 2) + std::pow(sigma(4, 0), 2) + std::pow(sigma(5, 0), 2)));
	//	//hexa.vms_vox(i) = vms;
	//}
}

void homo3D::FeSensi(hex_element& hexa)
{
	int node_size = hexa.periodic_cell_dofs.maxCoeff() + 1;
	int cell_size = hexa.cells_.rows();
	Eigen::VectorXd pdensity;
	pdensity.setZero(hexa.densities.rows());
	Eigen::VectorXd pdensity2;
	pdensity2.setZero(hexa.densities.rows());

	Eigen::MatrixXd Ke, Fe;
	getCellStiffnessAndForce(0, hexa.unit_mic_c, Ke, Fe, hexa);

	//Assembly
	Eigen::SparseMatrix<double> K(node_size - 3, node_size - 3);
	Eigen::MatrixXd F = Eigen::MatrixXd::Zero(node_size - 3, 6);
	{

		std::vector<Eigen::Triplet<double>> triplets;

		for (int i = 0; i < cell_size; ++i) {
			//Ke of element with different densities
			Eigen::MatrixXd cKe, cFe;
			pdensity(i) = pow(hexa.densities(i), hexa.panalty);
			pdensity2(i) = pow(hexa.densities(i), hexa.panalty - 1);

			cKe = Ke * (hexa.Emin + pdensity(i) * (hexa.Et - hexa.Emin));
			cFe = Fe * (hexa.Emin + pdensity(i) * (hexa.Et - hexa.Emin));

			for (int r = 0; r < 24; ++r) {
				int index_r = hexa.periodic_cell_dofs(i, r) - 3;
				if (index_r < 0) continue;

				for (int c = 0; c < 24; ++c) {
					int index_c = hexa.periodic_cell_dofs(i, c) - 3;
					if (index_c < 0) continue;
					if (Ke(r, c) != 0)
						triplets.emplace_back(index_r, index_c, cKe(r, c));
				}
				F.row(index_r) += cFe.row(r);
			}
		}

		K.setFromTriplets(triplets.begin(), triplets.end());
		K.makeCompressed();
	}

	//solve KX=f
	Eigen::MatrixXd X;
	{
		Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Upper> cg;
		//cg.setMaxIterations(200);
		cg.setTolerance(1e-6);
		cg.compute(K);
		X = cg.solve(F);
	}

	//fix dof= {0,1,2,4,5,11} for cell stiffness, however I do not know why...
	std::vector<int> unfixed_dofs = { 3,6,7,8,9,10,12,13,14,15,16,17,18,19,20,21,22,23 };

	Eigen::MatrixXd Ke_(unfixed_dofs.size(), unfixed_dofs.size()), Fe_(unfixed_dofs.size(), 6);

	for (int r = 0; r < unfixed_dofs.size(); ++r) {
		for (int c = 0; c < unfixed_dofs.size(); ++c) {
			Ke_(r, c) = Ke(unfixed_dofs[r], unfixed_dofs[c]);
		}
		Fe_.row(r) = Fe.row(unfixed_dofs[r]);
	}

	Eigen::MatrixXd temp = Ke_.lu().solve(Fe_);

	Eigen::MatrixXd X0 = Eigen::MatrixXd::Zero(24, 6);
	for (int i = 0; i < unfixed_dofs.size(); ++i)
		X0.row(unfixed_dofs[i]) = temp.row(i);


	//pre-computing
	std::vector<Eigen::MatrixXd> TX0(6), TX(6);
#pragma omp parallel for
	for (int i = 0; i < 6; i++) {
		TX0[i] = X0.col(i).transpose().replicate(cell_size, 1);

		TX[i].resize(cell_size, 24);

		for (int r = 0; r < cell_size; r++) {
			for (int c = 0; c < 24; c++) {
				TX[i](r, c) = hexa.periodic_cell_dofs(r, c) - 3 < 0 ? 0 : X(hexa.periodic_cell_dofs(r, c) - 3, i);
			}
		}
	}


	hexa.macro_constructive_tensor = Eigen::MatrixXd::Zero(6, 6);

	/*!< \f[C^H_{ij}=1/|\Omega|\sum_e\int_{\Omega_e}(\Chi_e^{0(i)}-\Chi_e^(i))^T k_e (\Chi_e^{0(j)}-\Chi_e^(j))d\Omega^(e)\f] */
	Eigen::MatrixXd left, right;
	Eigen::VectorXd n1;
	n1.setOnes(pdensity.rows());
	double volume = 1 / (hexa.nodes_.colwise().maxCoeff() - hexa.nodes_.colwise().minCoeff()).prod();

	hexa.dEH.clear();

	for (int i = 0; i < 6; i++) {
		left = (TX0[i] - TX[i]) * Ke;

		for (int j = 0; j < 6; j++) {
			right = TX0[j] - TX[j];

			Eigen::VectorXd sumrows = (left.array() * right.array()).rowwise().sum();
			hexa.macro_constructive_tensor(i, j) = sumrows.transpose() * (hexa.Emin * n1 + pdensity * (hexa.Et - hexa.Emin));
			n1.setOnes(pdensity.rows());
			//hexa.macro_constructive_tensor(i, j) = (left.array() * right.array()).sum();

			Eigen::VectorXd dEHij = hexa.panalty * (hexa.Et - hexa.Emin) * (pdensity2.array() * sumrows.array());

			dEHij *= volume;

			hexa.dEH.push_back(dEHij);
		}
	}

	hexa.macro_constructive_tensor *= volume;

	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 6; j++) {
			if (abs(hexa.macro_constructive_tensor(i, j)) < 1e-5)
			{
				hexa.macro_constructive_tensor(i, j) = 0;
			}
		}
	}
}

void homo3D::get_b(hex_element& hexa)
{
	hexa.B.resize(6, 24);
	double x = 0.0, y = 0.0, z = 0.0;
	Eigen::MatrixXd q(3, 8), corner(8, 3), jacobian(3, 3);
	Eigen::MatrixXd f_strain(6, 6);
	q << -((y - 1) * (z - 1)) / 8, ((y - 1) * (z - 1)) / 8,
		-((y + 1) * (z - 1)) / 8, ((y + 1) * (z - 1)) / 8,
		((y - 1) * (z + 1)) / 8, -((y - 1) * (z + 1)) / 8,
		((y + 1) * (z + 1)) / 8, -((y + 1) * (z + 1)) / 8,
		-((x - 1) * (z - 1)) / 8, ((x + 1) * (z - 1)) / 8,
		-((x + 1) * (z - 1)) / 8, ((x - 1) * (z - 1)) / 8,
		((x - 1) * (z + 1)) / 8, -((x + 1) * (z + 1)) / 8,
		((x + 1) * (z + 1)) / 8, -((x - 1) * (z + 1)) / 8,
		-((x - 1) * (y - 1)) / 8, ((x + 1) * (y - 1)) / 8,
		-((x + 1) * (y + 1)) / 8, ((x - 1) * (y + 1)) / 8,
		((x - 1) * (y - 1)) / 8, -((x + 1) * (y - 1)) / 8,
		((x + 1) * (y + 1)) / 8, -((x - 1) * (y + 1)) / 8;
	for (int m = 0; m < 8; ++m) {
		corner.row(m) = hexa.nodes_.row(hexa.cells_(0, m));
	}
	jacobian.noalias() = q * corner;

	Eigen::MatrixXd qxyz = jacobian.lu().solve(q);
	//elib::writeCSV("../data/qxyz.csv", qxyz, true);



	for (int m = 0; m < 8; m++) {
		Eigen::MatrixXd Be(6, 3);
		Be << qxyz(0, m), 0, 0,
			0, qxyz(1, m), 0,
			0, 0, qxyz(2, m),
			qxyz(1, m), qxyz(0, m), 0,
			0, qxyz(2, m), qxyz(1, m),
			qxyz(2, m), 0, qxyz(0, m);
		hexa.B.block(0, 0 + 3 * m, 6, 3) = Be;
	}
}
void homo3D::solve_lsm(hex_element& hexa, std::vector<Eigen::VectorXd>& velosity)
{
	int node_size = hexa.periodic_cell_dofs.maxCoeff() + 1;
	int cell_size = hexa.cells_.rows();

	Eigen::MatrixXd Ke, Fe;
	getCellStiffnessAndForce(0, hexa.micro_constructive_tensor, Ke, Fe, hexa);

	//Assembly
	Eigen::SparseMatrix<double> K(node_size - 3, node_size - 3);
	Eigen::MatrixXd F = Eigen::MatrixXd::Zero(node_size - 3, 6);
	{

		std::vector<Eigen::Triplet<double>> triplets;

		for (int i = 0; i < cell_size; ++i) {

			for (int r = 0; r < 24; ++r) {
				int index_r = hexa.periodic_cell_dofs(i, r) - 3;
				if (index_r < 0) continue;

				for (int c = 0; c < 24; ++c) {
					int index_c = hexa.periodic_cell_dofs(i, c) - 3;
					if (index_c < 0) continue;
					if (Ke(r, c) != 0)
						triplets.emplace_back(index_r, index_c, Ke(r, c));
				}
				F.row(index_r) += Fe.row(r);
			}
		}

		K.setFromTriplets(triplets.begin(), triplets.end());
		K.makeCompressed();
		//Eigen::MatrixXd K_ = K;
		//std::cout << K.block(0, 0, 50, 50) << std::endl;


		//std::cout << "Assemble time: " << end - start << " ms" << std::endl;
	}

	//solve KX=f
	auto start = std::chrono::high_resolution_clock::now();
	
	Eigen::MatrixXd X;
	{
		Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Upper> cg;
		//cg.setMaxIterations(200);
		cg.setTolerance(1e-6);
		cg.compute(K);
		X = cg.solve(F);
	}
	auto end = std::chrono::high_resolution_clock::now();
	std::chrono::duration<double> duration = end - start;
	std::cout << "KX=f运行时间: " << duration.count() << " 秒." << std::endl;

	//fix dof= {0,1,2,4,5,11} for cell stiffness, however I do not know why...
	std::vector<int> unfixed_dofs = { 3,6,7,8,9,10,12,13,14,15,16,17,18,19,20,21,22,23 };

	Eigen::MatrixXd Ke_(unfixed_dofs.size(), unfixed_dofs.size()), Fe_(unfixed_dofs.size(), 6);

	for (int r = 0; r < unfixed_dofs.size(); ++r) {
		for (int c = 0; c < unfixed_dofs.size(); ++c) {
			Ke_(r, c) = Ke(unfixed_dofs[r], unfixed_dofs[c]);
		}
		Fe_.row(r) = Fe.row(unfixed_dofs[r]);
	}

	Eigen::MatrixXd temp = Ke_.lu().solve(Fe_);

	Eigen::MatrixXd X0 = Eigen::MatrixXd::Zero(24, 6);
	for (int i = 0; i < unfixed_dofs.size(); ++i)
		X0.row(unfixed_dofs[i]) = temp.row(i);

	//pre-computing
	std::vector<Eigen::MatrixXd> TX0(6), TX(6);
	for (int i = 0; i < 6; i++) {
		TX0[i] = X0.col(i).transpose().replicate(cell_size, 1);

		TX[i].resize(cell_size, 24);

		for (int r = 0; r < cell_size; r++) {
			for (int c = 0; c < 24; c++) {
				TX[i](r, c) = hexa.periodic_cell_dofs(r, c) - 3 < 0 ? 0 : X(hexa.periodic_cell_dofs(r, c) - 3, i);
			}
		}
	}

	hexa.macro_constructive_tensor = Eigen::MatrixXd::Zero(6, 6);

	double volume = (hexa.nodes_.colwise().maxCoeff() - hexa.nodes_.colwise().minCoeff()).prod();
	//std::cout << "volume: " << volume << std::endl;

	for (int i = 0; i < velosity.size(); i++) {
		Eigen::MatrixXd dcht = Eigen::MatrixXd::Zero(6, 6);
		hexa.dCH.push_back(dcht);	
	}

	/*!< \f[C^H_{ij}=1/|\Omega|\sum_e\int_{\Omega_e}(\Chi_e^{0(i)}-\Chi_e^(i))^T k_e (\Chi_e^{0(j)}-\Chi_e^(j))d\Omega^(e)\f] */
	Eigen::MatrixXd left, right;
	for (int i = 0; i < 6; i++) {
		left.noalias() = (TX0[i] - TX[i]) * Ke;

		for (int j = 0; j < 6; j++) {
			right.noalias() = TX0[j] - TX[j];
			Eigen::VectorXd sumrows = (left.array() * right.array()).rowwise().sum();
			

			hexa.macro_constructive_tensor(i, j) = (left.array() * right.array()).sum();

			if (abs(hexa.macro_constructive_tensor(i, j)) < 1e-5) hexa.macro_constructive_tensor(i, j) = 0;

			for (int k = 0; k < velosity.size(); k++) {
				hexa.dCH[k](i, j) = sumrows.transpose() * velosity[k];
				hexa.dCH[k](i, j) *= 1 / volume;
				if (abs(hexa.dCH[k](i, j)) < 1e-5) hexa.dCH[k](i, j) = 0;
			}

		}
	}

	hexa.macro_constructive_tensor *= 1 / (hexa.nodes_.colwise().maxCoeff() - hexa.nodes_.colwise().minCoeff()).prod();

	hexa.macro_compliance_tensor = hexa.macro_constructive_tensor.inverse();

	//多了一部分
	get_b(hexa);
	Eigen::MatrixXd d(24, 6);
	hexa.vmsForEnergyAbsorption = 0;
	for (int i = 0; i < cell_size; i++)
	{
		for (int j = 0; j < 24; j++)
		{
			for (int r = 0; r < 6; r++) {
				d(j, r) = hexa.periodic_cell_dofs(i, j) - 3 < 0 ? 0 : X(hexa.periodic_cell_dofs(i, j) - 3, r);
			}
		}

		Eigen::MatrixXd td = X0 - d;
		Eigen::MatrixXd sigma = hexa.micro_constructive_tensor * hexa.B * td;
		double vms = std::sqrt(0.5 * (std::pow((sigma(0, 0) - sigma(1, 0)), 2) + std::pow((sigma(1, 0) - sigma(2, 0)), 2) + std::pow((sigma(2, 0) - sigma(0, 0)), 2)) + 3 * (std::pow(sigma(3, 0), 2) + std::pow(sigma(4, 0), 2) + std::pow(sigma(5, 0), 2)));
		hexa.vmsForEnergyAbsorption += vms;
	}
	hexa.dVmsForEnergyAbsorption = 0;

}
