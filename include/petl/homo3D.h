#pragma once
#include <Eigen/Sparse>
#include <Eigen/Core>
#include <Eigen/Eigen>
#include <iostream>
#include <vector>
#include <ctime>
#include <set>
#include "write_csv.h"

namespace homo3D {

	//bool debug_out = false;
	struct hex_element {

		Eigen::MatrixXi periodic_cell_dofs;

		std::vector<Eigen::MatrixXd> dCH;//for lsm

		Eigen::MatrixXd micro_constructive_tensor;
		Eigen::MatrixXd macro_constructive_tensor;
		Eigen::MatrixXd macro_compliance_tensor;

		Eigen::MatrixXd unit_mic_c;
		Eigen::SparseMatrix<double> H;
		Eigen::VectorXd Hs;


		std::vector<Eigen::VectorXd> dEH;//for simp

		double Emin = 1e-9;
		//double Et = 1500;
		double Et = 1;
		double v = 0.3;
		
		double panalty = 3;
		double rmin = 1.5;

		Eigen::VectorXd densities;

		Eigen::MatrixXd nodes_;//coordinates
		Eigen::MatrixXi cells_;//index
		std::vector<Eigen::VectorXd> velosity;
		Eigen::VectorXd vms_vox;

		//比PSL多出来的
		double vmsForEnergyAbsorption;
		double dVmsForEnergyAbsorption;
		Eigen::MatrixXd B;
	};

	
	void init(const Eigen::MatrixXd& grids, const Eigen::VectorXi& labels, const Eigen::Vector3i& ranges, hex_element& hexa);

	void init_simp(const Eigen::MatrixXd& grids, const Eigen::VectorXi& labels, const Eigen::Vector3i& ranges, hex_element& hexa);

	void setMaterial(const double E, const double v, hex_element& hexa);

	void setMaterial_simp(const double v, hex_element& hexa);	//不管E了,默认为1

	void get_sigma_f(const int index, Eigen::VectorXd& sigma, Eigen::MatrixXd& f, Eigen::MatrixXd& d, Eigen::VectorXd& m1, hex_element& hexa);

	void getCellStiffnessAndForce(const int index, const Eigen::MatrixXd& construtive_tensor, Eigen::MatrixXd& stiffness, Eigen::MatrixXd& force, hex_element& hexa);

	//void getUnitCellStiffnessAndForce(const Eigen::MatrixXd& construtive_tensor, Eigen::MatrixXd& stiffness, Eigen::MatrixXd& force, hex_element& hexa);

	void FeSensi(hex_element& hexa);

	void solve(hex_element& hexa);

	void solve_lsm(hex_element& hexa, std::vector<Eigen::VectorXd>& velosity);

	void get_b(hex_element& hexa);
}