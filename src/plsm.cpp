#include "petl/plsm.h"
#include "petl/homo3D.h"
#include <MMASolver.h>
#include <cassert>
#include <igl/writeOBJ.h>
#include <filesystem>
#include <stdexcept>

void plsmtoSolid(Truss& PTL, double trd, int id, double begin_smooth_factor, const std::string& skeleton_path, int max_iterations)
{
	//默认参数
	double unit_size = PTL.cell_length;	//一个单胞的尺寸

	double relative_r_outer = 0.04;			//相对半径,然后再*unit_size
	double relative_smooth_factor = 0.01;

	double r_outer = unit_size * relative_r_outer; //范围1.5-4左右
	//double smooth_factor = unit_size * relative_smooth_factor;//(0,3)左右 不能取0
	double smooth_factor = begin_smooth_factor;

	Eigen::Vector3i range(PTL.resolution + 1, PTL.resolution + 1, PTL.resolution + 1);


	//读入输入文件并输出初始模型 需要修改
	PTL.setRInner(0, -1);
	PTL.setROuter(r_outer, -1);
	PTL.setSmoothFactor(smooth_factor, -1);
	//int id = 3;
	PTL.getSkeletonFromPanetta(skeleton_path);

	if (!PTL.setSolidVolumeFraction(trd))
		throw std::runtime_error("Unable to initialize the requested volume fraction");
	PTL.genSolidTrussVoxel();			//生成体素
	PTL.genSolidTrussMesh();				//生成模型
	PTL.genTrussGrid();				//生成homo3d用网格
	std::string folder_name = PTL.output_path + std::to_string(id) + "-trd" + std::to_string(trd) + "-maxEwithIso";
	std::filesystem::create_directories(folder_name);
	std::string obj_path = folder_name + "/" + std::to_string(id) + "-iter0.obj";
	std::string csv_path = folder_name + "/" + std::to_string(id) + "-iter0.csv";
	if (!igl::writeOBJ(obj_path, PTL.V, PTL.F))
		throw std::runtime_error("Cannot write mesh: " + obj_path);
	if (!elib::writeCSV(csv_path, PTL.voxel))
		throw std::runtime_error("Cannot write voxels: " + csv_path);

	std::cout << "初始体积分数为" << PTL.getVolumeFraction() << std::endl;


	//log文件
	std::string opt_path = folder_name + "/" + "opt.txt";
	std::ofstream file;
	file.open(opt_path);
	if (!file) throw std::runtime_error("Cannot write optimization log: " + opt_path);
	file << "-------------------- P L S M T O --------------------" << std::endl;

	//MMA初始化相关变量
	double change = 1;
	int num_iter = 0;

	//int n = 2;	//优化参数的个数 这里只有r_outer和smooth_factor
	//int m = 2;	//优化约束的个数
	int n = PTL.base_edges.size() + PTL.base_joints.size();
	int m = 3;

	double r_outer_min = 0.5;
	double r_outer_max = 6;
	double smooth_factor_min = 0.1;
	double smooth_factor_max = 4;

	double r_outer_movlim = 1.0;
	double smooth_factor_movlim = 0.3;

	std::vector<double> x_storage(n);
	double* x = x_storage.data();		//实际优化的所有参数组成变量
	std::vector<double> x_old_storage(n);
	double* x_old = x_old_storage.data();
	std::vector<double> f_storage(1);
	double* f = f_storage.data();		//目标量 如CH[1,1]
	std::vector<double> df_storage(n);
	double* df = df_storage.data();		//梯度
	std::vector<double> g_storage(m);
	double* g = g_storage.data();		//体积分数和各向同性约束
	std::vector<double> dg_storage(n * m);
	double* dg = dg_storage.data();	//体积分数的梯度和各向同性的梯度
	std::vector<double> x_min_storage(n);
	double* x_min = x_min_storage.data();
	std::vector<double> x_max_storage(n);
	double* x_max = x_max_storage.data();

	/*
	int x_index = 0;
	x[x_index] = PTL.r_outer;
	x_old[x_index] = PTL.r_outer;
	x_index++;

	x[x_index] = PTL.smooth_factor;
	x_old[x_index] = PTL.smooth_factor;
	x_index++;*/
	int x_index = 0;
	for (int i = 0; i < PTL.base_edges.size(); i++)
	{
		x[x_index] = PTL.base_edges[i].r_outer;
		x_old[x_index] = PTL.base_edges[i].r_outer;
		x_index++;
	}
	for (int i = 0; i < PTL.base_joints.size(); i++)
	{
		x[x_index] = PTL.base_joints[i].smooth_factor;
		x_old[x_index] = PTL.base_joints[i].smooth_factor;
		x_index++;
	}

	//MMA优化器初始化
	MMASolver mma(n, m);

	while (1)
	{
		homo3D::hex_element hex_e;
		setMaterial(1, 0.35, hex_e);
		init(PTL.grid, PTL.voxel, range, hex_e);
		std::vector<Eigen::VectorXd> velosity;

		PTL.getSolidVelocities(velosity);

		solve_lsm(hex_e, velosity);	//组合出最后梯度

		Eigen::MatrixXd CH = hex_e.macro_constructive_tensor;

		Eigen::MatrixXd SH = CH.inverse();
		//计算dEH
		Eigen::MatrixXd dEH;
		std::vector<Eigen::MatrixXd> dSHv;
		for (int i = 0; i < hex_e.dCH.size(); i++)
		{
			Eigen::MatrixXd dshi;
			dshi.setZero(6, 6);
			dSHv.push_back(dshi);
		}
		for (int i = 0; i < n; i++) {
			dEH.setZero(6, 6);
			for (int j = 0; j < 6; j++) {
				for (int k = 0; k < 6; k++) {
					dEH(j, k) = hex_e.dCH[i](j, k);
				}
			}
			Eigen::MatrixXd dSH = -1 * SH * dEH * SH;
			for (int j = 0; j < 6; j++) {
				for (int k = 0; k < 6; k++) {
					dSHv[i](j, k) = dSH(j, k);
				}
			}
		}
		//isotropy constraint
		double ehisoii = (CH(0, 0) + CH(1, 1) + CH(2, 2)) / 3;
		double ehisoij = (CH(0, 1) + CH(0, 2) + CH(1, 2)) / 3;
		double ehisoii2 = (ehisoii - ehisoij) / 2;
		Eigen::MatrixXd EH_iso(6, 6);
		EH_iso << ehisoii, ehisoij, ehisoij, 0, 0, 0,
			ehisoij, ehisoii, ehisoij, 0, 0, 0,
			ehisoij, ehisoij, ehisoii, 0, 0, 0,
			0, 0, 0, ehisoii2, 0, 0,
			0, 0, 0, 0, ehisoii2, 0,
			0, 0, 0, 0, 0, ehisoii2;
		for (int i = 0; i < 36; i++)
		{
			Eigen::VectorXd dshi;
			dshi.setZero(n);
			hex_e.dEH.push_back(dshi);
		}
		for (int j = 0; j < 6; j++) {
			for (int k = 0; k < 6; k++) {
				for (int i = 0; i < n; i++)
				{
					hex_e.dEH[j * 6 + k](i) = hex_e.dCH[i](j, k);
				}
			}
		}
		Eigen::VectorXd dehisoii = (hex_e.dEH[0] + hex_e.dEH[7] + hex_e.dEH[14]) / 3;
		Eigen::VectorXd dehisoij = (hex_e.dEH[1] + hex_e.dEH[2] + hex_e.dEH[8]) / 3;
		//Eigen::VectorXd dehisoii2 = (dehisoii + dehisoij) / 2;		//怎么变成+了 于2024.3.31修改
		Eigen::VectorXd dehisoii2 = (dehisoii - dehisoij) / 2;
		double g_iso = 0;
		Eigen::VectorXd dgiso = Eigen::VectorXd::Zero(n);
		for (int i = 0; i < 6; i++) {
			for (int j = 0; j < 6; j++) {
				g_iso += pow(EH_iso(i, j) - CH(i, j), 2);
				if (i < 3 && j < 3 && i == j)
				{
					dgiso += 2 * (EH_iso(i, j) - CH(i, j)) / EH_iso(0, 0) * ((dehisoii - hex_e.dEH[i * 6 + j]) * EH_iso(0, 0) - dehisoii * (EH_iso(i, j) - CH(i, j))) / (EH_iso(0, 0) * EH_iso(0, 0));
				}
				else if (i < 3 && j < 3 && i != j)
				{
					dgiso += 2 * (EH_iso(i, j) - CH(i, j)) / EH_iso(0, 0) * ((dehisoij - hex_e.dEH[i * 6 + j]) * EH_iso(0, 0) - dehisoii * (EH_iso(i, j) - CH(i, j))) / (EH_iso(0, 0) * EH_iso(0, 0));
				}
				else if (i > 2 && j > 2 && i == j)
				{
					dgiso += 2 * (EH_iso(i, j) - CH(i, j)) / EH_iso(0, 0) * ((dehisoii2 - hex_e.dEH[i * 6 + j]) * EH_iso(0, 0) - dehisoii * (EH_iso(i, j) - CH(i, j))) / (EH_iso(0, 0) * EH_iso(0, 0));
				}
				else
				{
					dgiso += 2 * (-CH(i, j)) / EH_iso(0, 0) * ((-hex_e.dEH[i * 6 + j]) * EH_iso(0, 0) - dehisoii * (-CH(i, j))) / (EH_iso(0, 0) * EH_iso(0, 0));
				}
			}
		}
		g_iso /= pow(EH_iso(0, 0), 2);

		//A_U计算
		double Kv = (CH(0, 0) + CH(1, 1) + CH(2, 2)) / 9.0 + (CH(0, 1) + CH(0, 2) + CH(1, 2)) * 2.0 / 9.0;
		double Gv = (CH(0, 0) + CH(1, 1) + CH(2, 2) - CH(0, 1) - CH(0, 2) - CH(1, 2)) / 15.0 + (CH(3, 3) + CH(4, 4) + CH(5, 5)) / 5.0;
		double Kr = 1.0 / (SH(0, 0) + SH(1, 1) + SH(2, 2) + 2.0 * SH(0, 1) + 2.0 * SH(0, 2) + 2.0 * SH(1, 2));
		double Gr = 15.0 / (4.0 * SH(0, 0) + 4.0 * SH(1, 1) + 4.0 * SH(2, 2) - 4.0 * SH(0, 1) - 4.0 * SH(0, 2) - 4.0 * SH(1, 2) + 3.0 * SH(3, 3) + 3.0 * SH(4, 4) + 3.0 * SH(5, 5));
		double Au = 5.0 * Gv / Gr + Kv / Kr - 6;


		double crd = (double)PTL.voxel.sum() / (double)PTL.voxel.rows();	//当前相对密度or体积分数

		//设置目标函数
		f[0] = -1 * (CH(0, 0) + CH(1, 1) + CH(2, 2));	//在尽可能各向同性的情况下max E
		//f[0] = CH(0, 1) + CH(0, 2) + CH(1, 2) - powf(0.5, num_iter) * (CH(0, 0) + CH(1, 1) + CH(2, 2)); //minNu
		//f[0] = -1 * (CH(3, 3) + CH(4, 4) + CH(5, 5));	//maxG
		//f[0] = Au;	//iso

		//f[0] = -1 * (CH(0, 0));						//max某方向上的E
		//f[0] = -1 * (CH(3, 3));						//max某方向上的G
		//f[0] = g_iso;


		std::cout << "iter_" << num_iter << " " << "rd: " << crd << " " << "obj: " << f[0] << " " << "g_iso: " << g_iso << std::endl;
		std::cout << CH << std::endl;

		std::cout << "E: " << 1.0 / SH(0, 0) << std::endl;
		std::cout << "nu: " << -1 * SH(1, 0) / SH(0, 0) << std::endl;
		std::cout << "G: " << 1.0 / SH(3, 3) << std::endl;
		std::cout << "A^U: " << Au << std::endl;
		std::cout << "Zener ratio: " << 2 * CH(3, 3) / (CH(0, 0) - CH(0, 1)) << std::endl;
		std::cout << std::endl;

		file << "iter_" << num_iter << " " << "rd: " << crd << " " << "obj: " << f[0] << " " << "g_iso: " << g_iso << std::endl << std::endl;
		file << CH << std::endl << std::endl;

		file << "E: " << 1.0 / SH(0, 0) << std::endl;
		file << "nu: " << -1 * SH(1, 0) / SH(0, 0) << std::endl;
		file << "G: " << 1.0 / SH(3, 3) << std::endl;
		file << "A^U: " << Au << std::endl;
		file << "Zener ratio: " << 2 * CH(3, 3) / (CH(0, 0) - CH(0, 1)) << std::endl;
		file << std::endl;

		if (!((num_iter < 10 || change > 1e-2) && num_iter < max_iterations))	//终止条件,最少迭代10次 最大迭代150次?
			break;

		if (crd > trd + 0.05)	//超出体积分数约束
		{
			std::cout << "Out of target volume fraction!" << std::endl;
			file << "Out of target volume fraction!" << std::endl;
			break;
		}

		//dCH.size()=n 即优化参数个数
		for (int i = 0; i < hex_e.dCH.size(); i++) {
			df[i] = -1 * (hex_e.dCH[i](0, 0) + hex_e.dCH[i](1, 1) + hex_e.dCH[i](2, 2));
			//df[i] = hex_e.dCH[i](0, 1) + hex_e.dCH[i](0, 2) + hex_e.dCH[i](1, 2) - powf(0.5, num_iter) * (hex_e.dCH[i](0, 0) + hex_e.dCH[i](1, 1) + hex_e.dCH[i](2, 2));
			//df[i] = -1 * (hex_e.dCH[i](3, 3) + hex_e.dCH[i](4, 4) + hex_e.dCH[i](5, 5));
			//df[i] = 1.0 / 45.0 * ((CH(0, 0) + CH(1, 1) + CH(2, 2) - CH(0, 1) - CH(0, 2) - CH(1, 2)) + 3 * (CH(3, 3) + CH(4, 4) + CH(5, 5))) * (4.0 * dSHv[i](0, 0) + 4.0 * dSHv[i](1, 1) + 4.0 * dSHv[i](2, 2) - 4.0 * dSHv[i](0, 1) - 4.0 * dSHv[i](0, 2) - 4.0 * dSHv[i](1, 2) + 3.0 * dSHv[i](3, 3) + 3.0 * dSHv[i](4, 4) + 3.0 * dSHv[i](5, 5))
			//	+ 1.0 / 45.0 * ((hex_e.dCH[i](0, 0) + hex_e.dCH[i](1, 1) + hex_e.dCH[i](2, 2) - hex_e.dCH[i](0, 1) - hex_e.dCH[i](0, 2) - hex_e.dCH[i](1, 2)) + 3 * (hex_e.dCH[i](3, 3) + hex_e.dCH[i](4, 4) + hex_e.dCH[i](5, 5))) * (4.0 * SH(0, 0) + 4.0 * SH(1, 1) + 4.0 * SH(2, 2) - 4.0 * SH(0, 1) - 4.0 * SH(0, 2) - 4.0 * SH(1, 2) + 3.0 * SH(3, 3) + 3.0 * SH(4, 4) + 3.0 * SH(5, 5))
			//	+ 1.0 / 9.0 * (CH(0, 0) + CH(1, 1) + CH(2, 2) + 2.0 * (CH(0, 1) + CH(0, 2) + CH(1, 2))) * (dSHv[i](0, 0) + dSHv[i](1, 1) + dSHv[i](2, 2) + 2.0 * dSHv[i](0, 1) + 2.0 * dSHv[i](0, 2) + 2.0 * dSHv[i](1, 2))
			//	+ 1.0 / 9.0 * (hex_e.dCH[i](0, 0) + hex_e.dCH[i](1, 1) + hex_e.dCH[i](2, 2) + 2.0 * (hex_e.dCH[i](0, 1) + hex_e.dCH[i](0, 2) + hex_e.dCH[i](1, 2))) * (SH(0, 0) + SH(1, 1) + SH(2, 2) + 2.0 * SH(0, 1) + 2.0 * SH(0, 2) + 2.0 * SH(1, 2));

			//df[i] = -1 * (hex_e.dCH[i](0, 0));
			//df[i] = -1 * (hex_e.dCH[i](3, 3));
			//df[i] = dgiso[i];
			std::cout << "df[" << i << "]=" << df[i] << std::endl;

			dg[i * m] = PTL.dV[i] / trd;		//体积分数偏导
			std::cout << "dg[" << i * m << "]=" << dg[i * m] << std::endl;

			dg[i * m + 1] = -1 * PTL.dV[i] / trd;	//体积分数偏导
			std::cout << "dg[" << i * m + 1 << "]=" << dg[i * m + 1] << std::endl;

			//dg[i * m + 2] = dgiso[i];	//iso偏导 PSL版
			//std::cout << "dg[" << i * m + 2 << "]=" << dg[i * m + 2] << std::endl;

			dg[i * m + 2] = 1.0 / 45.0 * ((CH(0, 0) + CH(1, 1) + CH(2, 2) - CH(0, 1) - CH(0, 2) - CH(1, 2)) + 3 * (CH(3, 3) + CH(4, 4) + CH(5, 5))) * (4.0 * dSHv[i](0, 0) + 4.0 * dSHv[i](1, 1) + 4.0 * dSHv[i](2, 2) - 4.0 * dSHv[i](0, 1) - 4.0 * dSHv[i](0, 2) - 4.0 * dSHv[i](1, 2) + 3.0 * dSHv[i](3, 3) + 3.0 * dSHv[i](4, 4) + 3.0 * dSHv[i](5, 5))
				+ 1.0 / 45.0 * ((hex_e.dCH[i](0, 0) + hex_e.dCH[i](1, 1) + hex_e.dCH[i](2, 2) - hex_e.dCH[i](0, 1) - hex_e.dCH[i](0, 2) - hex_e.dCH[i](1, 2)) + 3 * (hex_e.dCH[i](3, 3) + hex_e.dCH[i](4, 4) + hex_e.dCH[i](5, 5))) * (4.0 * SH(0, 0) + 4.0 * SH(1, 1) + 4.0 * SH(2, 2) - 4.0 * SH(0, 1) - 4.0 * SH(0, 2) - 4.0 * SH(1, 2) + 3.0 * SH(3, 3) + 3.0 * SH(4, 4) + 3.0 * SH(5, 5))
				+ 1.0 / 9.0 * (CH(0, 0) + CH(1, 1) + CH(2, 2) + 2.0 * (CH(0, 1) + CH(0, 2) + CH(1, 2))) * (dSHv[i](0, 0) + dSHv[i](1, 1) + dSHv[i](2, 2) + 2.0 * dSHv[i](0, 1) + 2.0 * dSHv[i](0, 2) + 2.0 * dSHv[i](1, 2))
				+ 1.0 / 9.0 * (hex_e.dCH[i](0, 0) + hex_e.dCH[i](1, 1) + hex_e.dCH[i](2, 2) + 2.0 * (hex_e.dCH[i](0, 1) + hex_e.dCH[i](0, 2) + hex_e.dCH[i](1, 2))) * (SH(0, 0) + SH(1, 1) + SH(2, 2) + 2.0 * SH(0, 1) + 2.0 * SH(0, 2) + 2.0 * SH(1, 2));

		}
		g[0] = crd / trd - 1;
		g[1] = 0.98 - crd / trd;
		//g[2] = g_iso - 1e-5;
		g[2] = Au - 1e-5;
		//std::cout << "g[0]=" << g[0] << ",g[1]=" << g[1] << std::endl;
		std::cout << "g[0]=" << g[0] << ",g[1]=" << g[1] << ",g[2]=" << g[2] << std::endl;

		//设置每次步进的限制范围
		x_index = 0;
		for (int i = 0; i < PTL.base_edges.size(); i++)
		{
			x_max[x_index] = std::min(r_outer_max, x[x_index] + r_outer_movlim);
			x_min[x_index] = std::max(r_outer_min, x[x_index] - r_outer_movlim);
			x_index++;
		}
		for (int i = 0; i < PTL.base_joints.size(); i++)
		{
			x_max[x_index] = std::min(smooth_factor_max, x[x_index] + smooth_factor_movlim);
			x_min[x_index] = std::max(smooth_factor_min, x[x_index] - smooth_factor_movlim);
			x_index++;
		}


		mma.Update(x, df, g, dg, x_min, x_max);		//求解



		//更新参数
		x_index = 0;
		for (int i = 0; i < PTL.base_edges.size(); i++)
		{
			std::cout << "R_outer[" << i << "]:" << x_old[x_index] << " -> " << x[x_index] << std::endl;
			file << "R_outer[" << i << "]:" << x_old[x_index] << " -> " << x[x_index] << std::endl;
			PTL.setROuter(x[x_index], i);
			x_index++;
		}
		for (int i = 0; i < PTL.base_joints.size(); i++)
		{
			std::cout << "smooth_factor[" << i << "]:" << x_old[x_index] << " -> " << x[x_index] << std::endl;
			file << "smooth_factor[" << i << "]:" << x_old[x_index] << " -> " << x[x_index] << std::endl;
			PTL.setSmoothFactor(x[x_index], i);
			x_index++;
		}

		double maxdiff = 1e-10;
		for (int i = 0; i < n; i++) {
			double diff = std::abs(x[i] - x_old[i]);
			maxdiff = std::max(maxdiff, diff);
			x_old[i] = x[i];
		}
		change = maxdiff;

		std::cout << "change:" << change << std::endl;
		std::cout << std::endl;
		file << "change:" << change << std::endl;
		file << std::endl;

		//重新生成
		num_iter++;
		PTL.genSolidTrussVoxel();			//生成体素
		PTL.genSolidTrussMesh();				//生成模型
		PTL.genTrussGrid();				//生成homo3d用网格
		obj_path = folder_name + "/" + std::to_string(id) + "-iter" + std::to_string(num_iter) + ".obj";
		csv_path = folder_name + "/" + std::to_string(id) + "-iter" + std::to_string(num_iter) + ".csv";
		if (!igl::writeOBJ(obj_path, PTL.V, PTL.F))
			throw std::runtime_error("Cannot write mesh: " + obj_path);
		if (!elib::writeCSV(csv_path, PTL.voxel))
			throw std::runtime_error("Cannot write voxels: " + csv_path);
	}

	file.close();
}

//注意outer_sdf和inner_sdf即可
void plsmtoHollow(Truss& PTL, double trd)
{

}

void plsmtoSolidETH(Truss& PTL, double trd, int id, double begin_smooth_factor, const std::string& skeleton_path, int max_iterations)
{
	//默认参数
	double unit_size = PTL.cell_length;	//一个单胞的尺寸

	double relative_r_outer = 0.04;			//相对半径,然后再*unit_size
	double relative_smooth_factor = 0.01;

	double r_outer = unit_size * relative_r_outer; //范围1.5-4左右
	//double smooth_factor = unit_size * relative_smooth_factor;//(0,3)左右 不能取0
	double smooth_factor = begin_smooth_factor;

	Eigen::Vector3i range(PTL.resolution + 1, PTL.resolution + 1, PTL.resolution + 1);


	//读入输入文件并输出初始模型 需要修改
	PTL.setRInner(0, -1);
	PTL.setROuter(r_outer, -1);
	PTL.setSmoothFactor(smooth_factor, -1);
	//int id = 3;
	PTL.getSkeletonFromETH(skeleton_path);

	if (!PTL.setSolidVolumeFractionETH(trd))
		throw std::runtime_error("Unable to initialize the requested volume fraction");
	PTL.genSolidTrussVoxelETH();			//生成体素
	PTL.genSolidTrussMeshETH();				//生成模型
	PTL.genTrussGrid();				//生成homo3d用网格
	std::string folder_name = PTL.output_path + std::to_string(id) + "-ETHtrd" + std::to_string(trd) + "-maxE";
	std::filesystem::create_directories(folder_name);
	std::string obj_path = folder_name + "/" + std::to_string(id) + "-iter0.obj";
	std::string csv_path = folder_name + "/" + std::to_string(id) + "-iter0.csv";
	if (!igl::writeOBJ(obj_path, PTL.V, PTL.F))
		throw std::runtime_error("Cannot write mesh: " + obj_path);
	if (!elib::writeCSV(csv_path, PTL.voxel))
		throw std::runtime_error("Cannot write voxels: " + csv_path);

	std::cout << "初始体积分数为" << PTL.getVolumeFraction() << std::endl;


	//log文件
	std::string opt_path = folder_name + "/" + "opt.txt";
	std::ofstream file;
	file.open(opt_path);
	if (!file) throw std::runtime_error("Cannot write optimization log: " + opt_path);
	file << "-------------------- P L S M T O --------------------" << std::endl;

	//MMA初始化相关变量
	double change = 1;
	int num_iter = 0;

	//int n = 2;	//优化参数的个数 这里只有r_outer和smooth_factor
	//int m = 2;	//优化约束的个数
	int n = PTL.base_edges.size() + PTL.base_joints.size();
	int m = 2;

	double r_outer_min = 0.5;
	double r_outer_max = 6;
	double smooth_factor_min = 0.1;
	double smooth_factor_max = 4;

	double r_outer_movlim = 1.0;
	double smooth_factor_movlim = 0.3;

	std::vector<double> x_storage(n);
	double* x = x_storage.data();		//实际优化的所有参数组成变量
	std::vector<double> x_old_storage(n);
	double* x_old = x_old_storage.data();
	std::vector<double> f_storage(1);
	double* f = f_storage.data();		//目标量 如CH[1,1]
	std::vector<double> df_storage(n);
	double* df = df_storage.data();		//梯度
	std::vector<double> g_storage(m);
	double* g = g_storage.data();		//体积分数和各向同性约束
	std::vector<double> dg_storage(n * m);
	double* dg = dg_storage.data();	//体积分数的梯度和各向同性的梯度
	std::vector<double> x_min_storage(n);
	double* x_min = x_min_storage.data();
	std::vector<double> x_max_storage(n);
	double* x_max = x_max_storage.data();

	/*
	int x_index = 0;
	x[x_index] = PTL.r_outer;
	x_old[x_index] = PTL.r_outer;
	x_index++;

	x[x_index] = PTL.smooth_factor;
	x_old[x_index] = PTL.smooth_factor;
	x_index++;*/
	int x_index = 0;
	for (int i = 0; i < PTL.base_edges.size(); i++)
	{
		x[x_index] = PTL.base_edges[i].r_outer;
		x_old[x_index] = PTL.base_edges[i].r_outer;
		x_index++;
	}
	for (int i = 0; i < PTL.base_joints.size(); i++)
	{
		x[x_index] = PTL.base_joints[i].smooth_factor;
		x_old[x_index] = PTL.base_joints[i].smooth_factor;
		x_index++;
	}

	//MMA优化器初始化
	MMASolver mma(n, m);

	while (1)
	{
		homo3D::hex_element hex_e;
		setMaterial(1, 0.35, hex_e);
		init(PTL.grid, PTL.voxel, range, hex_e);
		std::vector<Eigen::VectorXd> velosity;

		PTL.getSolidVelocitiesETH(velosity);

		solve_lsm(hex_e, velosity);	//组合出最后梯度

		Eigen::MatrixXd CH = hex_e.macro_constructive_tensor;

		Eigen::MatrixXd SH = CH.inverse();
		//计算dEH
		Eigen::MatrixXd dEH;
		std::vector<Eigen::MatrixXd> dSHv;
		for (int i = 0; i < hex_e.dCH.size(); i++)
		{
			Eigen::MatrixXd dshi;
			dshi.setZero(6, 6);
			dSHv.push_back(dshi);
		}
		for (int i = 0; i < n; i++) {
			dEH.setZero(6, 6);
			for (int j = 0; j < 6; j++) {
				for (int k = 0; k < 6; k++) {
					dEH(j, k) = hex_e.dCH[i](j, k);
				}
			}
			Eigen::MatrixXd dSH = -1 * SH * dEH * SH;
			for (int j = 0; j < 6; j++) {
				for (int k = 0; k < 6; k++) {
					dSHv[i](j, k) = dSH(j, k);
				}
			}
		}
		//isotropy constraint
		double ehisoii = (CH(0, 0) + CH(1, 1) + CH(2, 2)) / 3;
		double ehisoij = (CH(0, 1) + CH(0, 2) + CH(1, 2)) / 3;
		double ehisoii2 = (ehisoii - ehisoij) / 2;
		Eigen::MatrixXd EH_iso(6, 6);
		EH_iso << ehisoii, ehisoij, ehisoij, 0, 0, 0,
			ehisoij, ehisoii, ehisoij, 0, 0, 0,
			ehisoij, ehisoij, ehisoii, 0, 0, 0,
			0, 0, 0, ehisoii2, 0, 0,
			0, 0, 0, 0, ehisoii2, 0,
			0, 0, 0, 0, 0, ehisoii2;
		for (int i = 0; i < 36; i++)
		{
			Eigen::VectorXd dshi;
			dshi.setZero(n);
			hex_e.dEH.push_back(dshi);
		}
		for (int j = 0; j < 6; j++) {
			for (int k = 0; k < 6; k++) {
				for (int i = 0; i < n; i++)
				{
					hex_e.dEH[j * 6 + k](i) = hex_e.dCH[i](j, k);
				}
			}
		}
		Eigen::VectorXd dehisoii = (hex_e.dEH[0] + hex_e.dEH[7] + hex_e.dEH[14]) / 3;
		Eigen::VectorXd dehisoij = (hex_e.dEH[1] + hex_e.dEH[2] + hex_e.dEH[8]) / 3;
		//Eigen::VectorXd dehisoii2 = (dehisoii + dehisoij) / 2;		//怎么变成+了 于2024.3.31修改
		Eigen::VectorXd dehisoii2 = (dehisoii - dehisoij) / 2;
		double g_iso = 0;
		Eigen::VectorXd dgiso = Eigen::VectorXd::Zero(n);
		for (int i = 0; i < 6; i++) {
			for (int j = 0; j < 6; j++) {
				g_iso += pow(EH_iso(i, j) - CH(i, j), 2);
				if (i < 3 && j < 3 && i == j)
				{
					dgiso += 2 * (EH_iso(i, j) - CH(i, j)) / EH_iso(0, 0) * ((dehisoii - hex_e.dEH[i * 6 + j]) * EH_iso(0, 0) - dehisoii * (EH_iso(i, j) - CH(i, j))) / (EH_iso(0, 0) * EH_iso(0, 0));
				}
				else if (i < 3 && j < 3 && i != j)
				{
					dgiso += 2 * (EH_iso(i, j) - CH(i, j)) / EH_iso(0, 0) * ((dehisoij - hex_e.dEH[i * 6 + j]) * EH_iso(0, 0) - dehisoii * (EH_iso(i, j) - CH(i, j))) / (EH_iso(0, 0) * EH_iso(0, 0));
				}
				else if (i > 2 && j > 2 && i == j)
				{
					dgiso += 2 * (EH_iso(i, j) - CH(i, j)) / EH_iso(0, 0) * ((dehisoii2 - hex_e.dEH[i * 6 + j]) * EH_iso(0, 0) - dehisoii * (EH_iso(i, j) - CH(i, j))) / (EH_iso(0, 0) * EH_iso(0, 0));
				}
				else
				{
					dgiso += 2 * (-CH(i, j)) / EH_iso(0, 0) * ((-hex_e.dEH[i * 6 + j]) * EH_iso(0, 0) - dehisoii * (-CH(i, j))) / (EH_iso(0, 0) * EH_iso(0, 0));
				}
			}
		}
		g_iso /= pow(EH_iso(0, 0), 2);

		//A_U计算
		double Kv = (CH(0, 0) + CH(1, 1) + CH(2, 2)) / 9.0 + (CH(0, 1) + CH(0, 2) + CH(1, 2)) * 2.0 / 9.0;
		double Gv = (CH(0, 0) + CH(1, 1) + CH(2, 2) - CH(0, 1) - CH(0, 2) - CH(1, 2)) / 15.0 + (CH(3, 3) + CH(4, 4) + CH(5, 5)) / 5.0;
		double Kr = 1.0 / (SH(0, 0) + SH(1, 1) + SH(2, 2) + 2.0 * SH(0, 1) + 2.0 * SH(0, 2) + 2.0 * SH(1, 2));
		double Gr = 15.0 / (4.0 * SH(0, 0) + 4.0 * SH(1, 1) + 4.0 * SH(2, 2) - 4.0 * SH(0, 1) - 4.0 * SH(0, 2) - 4.0 * SH(1, 2) + 3.0 * SH(3, 3) + 3.0 * SH(4, 4) + 3.0 * SH(5, 5));
		double Au = 5.0 * Gv / Gr + Kv / Kr - 6;


		double crd = (double)PTL.voxel.sum() / (double)PTL.voxel.rows();	//当前相对密度or体积分数

		//设置目标函数
		f[0] = -1 * (CH(0, 0) + CH(1, 1) + CH(2, 2));	//在尽可能各向同性的情况下max E
		//f[0] = CH(0, 1) + CH(0, 2) + CH(1, 2) - powf(0.5, num_iter) * (CH(0, 0) + CH(1, 1) + CH(2, 2)); //minNu
		//f[0] = -1 * (CH(3, 3) + CH(4, 4) + CH(5, 5));	//maxG
		//f[0] = Au;	//iso

		//f[0] = -1 * (CH(0, 0));						//max某方向上的E
		//f[0] = -1 * (CH(3, 3));						//max某方向上的G
		//f[0] = g_iso;


		std::cout << "iter_" << num_iter << " " << "rd: " << crd << " " << "obj: " << f[0] << " " << "g_iso: " << g_iso << std::endl;
		std::cout << CH << std::endl;

		std::cout << "E: " << 1.0 / SH(0, 0) << std::endl;
		std::cout << "nu: " << -1 * SH(1, 0) / SH(0, 0) << std::endl;
		std::cout << "G: " << 1.0 / SH(3, 3) << std::endl;
		std::cout << "A^U: " << Au << std::endl;
		std::cout << "Zener ratio: " << 2 * CH(3, 3) / (CH(0, 0) - CH(0, 1)) << std::endl;
		std::cout << std::endl;

		file << "iter_" << num_iter << " " << "rd: " << crd << " " << "obj: " << f[0] << " " << "g_iso: " << g_iso << std::endl << std::endl;
		file << CH << std::endl << std::endl;

		file << "E: " << 1.0 / SH(0, 0) << std::endl;
		file << "nu: " << -1 * SH(1, 0) / SH(0, 0) << std::endl;
		file << "G: " << 1.0 / SH(3, 3) << std::endl;
		file << "A^U: " << Au << std::endl;
		file << "Zener ratio: " << 2 * CH(3, 3) / (CH(0, 0) - CH(0, 1)) << std::endl;
		file << std::endl;

		if (!((num_iter < 10 || change > 1e-2) && num_iter < max_iterations))	//终止条件,最少迭代10次 最大迭代150次?
			break;

		if (crd > trd + 0.05)	//超出体积分数约束
		{
			std::cout << "Out of target volume fraction!" << std::endl;
			file << "Out of target volume fraction!" << std::endl;
			break;
		}

		//dCH.size()=n 即优化参数个数
		for (int i = 0; i < hex_e.dCH.size(); i++) {
			df[i] = -1 * (hex_e.dCH[i](0, 0) + hex_e.dCH[i](1, 1) + hex_e.dCH[i](2, 2));
			//df[i] = hex_e.dCH[i](0, 1) + hex_e.dCH[i](0, 2) + hex_e.dCH[i](1, 2) - powf(0.5, num_iter) * (hex_e.dCH[i](0, 0) + hex_e.dCH[i](1, 1) + hex_e.dCH[i](2, 2));
			//df[i] = -1 * (hex_e.dCH[i](3, 3) + hex_e.dCH[i](4, 4) + hex_e.dCH[i](5, 5));
			//df[i] = 1.0 / 45.0 * ((CH(0, 0) + CH(1, 1) + CH(2, 2) - CH(0, 1) - CH(0, 2) - CH(1, 2)) + 3 * (CH(3, 3) + CH(4, 4) + CH(5, 5))) * (4.0 * dSHv[i](0, 0) + 4.0 * dSHv[i](1, 1) + 4.0 * dSHv[i](2, 2) - 4.0 * dSHv[i](0, 1) - 4.0 * dSHv[i](0, 2) - 4.0 * dSHv[i](1, 2) + 3.0 * dSHv[i](3, 3) + 3.0 * dSHv[i](4, 4) + 3.0 * dSHv[i](5, 5))
			//	+ 1.0 / 45.0 * ((hex_e.dCH[i](0, 0) + hex_e.dCH[i](1, 1) + hex_e.dCH[i](2, 2) - hex_e.dCH[i](0, 1) - hex_e.dCH[i](0, 2) - hex_e.dCH[i](1, 2)) + 3 * (hex_e.dCH[i](3, 3) + hex_e.dCH[i](4, 4) + hex_e.dCH[i](5, 5))) * (4.0 * SH(0, 0) + 4.0 * SH(1, 1) + 4.0 * SH(2, 2) - 4.0 * SH(0, 1) - 4.0 * SH(0, 2) - 4.0 * SH(1, 2) + 3.0 * SH(3, 3) + 3.0 * SH(4, 4) + 3.0 * SH(5, 5))
			//	+ 1.0 / 9.0 * (CH(0, 0) + CH(1, 1) + CH(2, 2) + 2.0 * (CH(0, 1) + CH(0, 2) + CH(1, 2))) * (dSHv[i](0, 0) + dSHv[i](1, 1) + dSHv[i](2, 2) + 2.0 * dSHv[i](0, 1) + 2.0 * dSHv[i](0, 2) + 2.0 * dSHv[i](1, 2))
			//	+ 1.0 / 9.0 * (hex_e.dCH[i](0, 0) + hex_e.dCH[i](1, 1) + hex_e.dCH[i](2, 2) + 2.0 * (hex_e.dCH[i](0, 1) + hex_e.dCH[i](0, 2) + hex_e.dCH[i](1, 2))) * (SH(0, 0) + SH(1, 1) + SH(2, 2) + 2.0 * SH(0, 1) + 2.0 * SH(0, 2) + 2.0 * SH(1, 2));

			//df[i] = -1 * (hex_e.dCH[i](0, 0));
			//df[i] = -1 * (hex_e.dCH[i](3, 3));
			//df[i] = dgiso[i];
			std::cout << "df[" << i << "]=" << df[i] << std::endl;

			dg[i * m] = PTL.dV[i] / trd;		//体积分数偏导
			std::cout << "dg[" << i * m << "]=" << dg[i * m] << std::endl;

			dg[i * m + 1] = -1 * PTL.dV[i] / trd;	//体积分数偏导
			std::cout << "dg[" << i * m + 1 << "]=" << dg[i * m + 1] << std::endl;

			//dg[i * m + 2] = dgiso[i];	//iso偏导 PSL版
			//std::cout << "dg[" << i * m + 2 << "]=" << dg[i * m + 2] << std::endl;

			/*dg[i * m + 2] = 1.0 / 45.0 * ((CH(0, 0) + CH(1, 1) + CH(2, 2) - CH(0, 1) - CH(0, 2) - CH(1, 2)) + 3 * (CH(3, 3) + CH(4, 4) + CH(5, 5))) * (4.0 * dSHv[i](0, 0) + 4.0 * dSHv[i](1, 1) + 4.0 * dSHv[i](2, 2) - 4.0 * dSHv[i](0, 1) - 4.0 * dSHv[i](0, 2) - 4.0 * dSHv[i](1, 2) + 3.0 * dSHv[i](3, 3) + 3.0 * dSHv[i](4, 4) + 3.0 * dSHv[i](5, 5))
				+ 1.0 / 45.0 * ((hex_e.dCH[i](0, 0) + hex_e.dCH[i](1, 1) + hex_e.dCH[i](2, 2) - hex_e.dCH[i](0, 1) - hex_e.dCH[i](0, 2) - hex_e.dCH[i](1, 2)) + 3 * (hex_e.dCH[i](3, 3) + hex_e.dCH[i](4, 4) + hex_e.dCH[i](5, 5))) * (4.0 * SH(0, 0) + 4.0 * SH(1, 1) + 4.0 * SH(2, 2) - 4.0 * SH(0, 1) - 4.0 * SH(0, 2) - 4.0 * SH(1, 2) + 3.0 * SH(3, 3) + 3.0 * SH(4, 4) + 3.0 * SH(5, 5))
				+ 1.0 / 9.0 * (CH(0, 0) + CH(1, 1) + CH(2, 2) + 2.0 * (CH(0, 1) + CH(0, 2) + CH(1, 2))) * (dSHv[i](0, 0) + dSHv[i](1, 1) + dSHv[i](2, 2) + 2.0 * dSHv[i](0, 1) + 2.0 * dSHv[i](0, 2) + 2.0 * dSHv[i](1, 2))
				+ 1.0 / 9.0 * (hex_e.dCH[i](0, 0) + hex_e.dCH[i](1, 1) + hex_e.dCH[i](2, 2) + 2.0 * (hex_e.dCH[i](0, 1) + hex_e.dCH[i](0, 2) + hex_e.dCH[i](1, 2))) * (SH(0, 0) + SH(1, 1) + SH(2, 2) + 2.0 * SH(0, 1) + 2.0 * SH(0, 2) + 2.0 * SH(1, 2));
			*/
		}
		g[0] = crd / trd - 1;
		g[1] = 0.98 - crd / trd;
		//g[2] = g_iso - 1e-5;
		//g[2] = Au - 1e-5;
		std::cout << "g[0]=" << g[0] << ",g[1]=" << g[1] << std::endl;
		//std::cout << "g[0]=" << g[0] << ",g[1]=" << g[1] << ",g[2]=" << g[2] << std::endl;

		//设置每次步进的限制范围
		x_index = 0;
		for (int i = 0; i < PTL.base_edges.size(); i++)
		{
			x_max[x_index] = std::min(r_outer_max, x[x_index] + r_outer_movlim);
			x_min[x_index] = std::max(r_outer_min, x[x_index] - r_outer_movlim);
			x_index++;
		}
		for (int i = 0; i < PTL.base_joints.size(); i++)
		{
			x_max[x_index] = std::min(smooth_factor_max, x[x_index] + smooth_factor_movlim);
			x_min[x_index] = std::max(smooth_factor_min, x[x_index] - smooth_factor_movlim);
			x_index++;
		}


		mma.Update(x, df, g, dg, x_min, x_max);		//求解



		//更新参数
		x_index = 0;
		for (int i = 0; i < PTL.base_edges.size(); i++)
		{
			std::cout << "R_outer[" << i << "]:" << x_old[x_index] << " -> " << x[x_index] << std::endl;
			file << "R_outer[" << i << "]:" << x_old[x_index] << " -> " << x[x_index] << std::endl;
			PTL.setROuter(x[x_index], i);
			x_index++;
		}
		for (int i = 0; i < PTL.base_joints.size(); i++)
		{
			std::cout << "smooth_factor[" << i << "]:" << x_old[x_index] << " -> " << x[x_index] << std::endl;
			file << "smooth_factor[" << i << "]:" << x_old[x_index] << " -> " << x[x_index] << std::endl;
			PTL.setSmoothFactor(x[x_index], i);
			x_index++;
		}

		double maxdiff = 1e-10;
		for (int i = 0; i < n; i++) {
			double diff = std::abs(x[i] - x_old[i]);
			maxdiff = std::max(maxdiff, diff);
			x_old[i] = x[i];
		}
		change = maxdiff;

		std::cout << "change:" << change << std::endl;
		std::cout << std::endl;
		file << "change:" << change << std::endl;
		file << std::endl;

		//重新生成
		num_iter++;
		PTL.genSolidTrussVoxelETH();			//生成体素
		PTL.genSolidTrussMeshETH();				//生成模型
		PTL.genTrussGrid();				//生成homo3d用网格
		obj_path = folder_name + "/" + std::to_string(id) + "-iter" + std::to_string(num_iter) + ".obj";
		csv_path = folder_name + "/" + std::to_string(id) + "-iter" + std::to_string(num_iter) + ".csv";
		if (!igl::writeOBJ(obj_path, PTL.V, PTL.F))
			throw std::runtime_error("Cannot write mesh: " + obj_path);
		if (!elib::writeCSV(csv_path, PTL.voxel))
			throw std::runtime_error("Cannot write voxels: " + csv_path);
	}

	file.close();
}

void plsmtoSolidModified(Truss& PTL, double trd, int id, double begin_smooth_factor, const std::string& skeleton_path, int max_iterations)
{
	//默认参数
	double unit_size = PTL.cell_length;	//一个单胞的尺寸

	double relative_r_outer = 0.04;			//相对半径,然后再*unit_size
	double relative_smooth_factor = 0.01;

	double r_outer = unit_size * relative_r_outer; //范围1.5-4左右
	//double smooth_factor = unit_size * relative_smooth_factor;//(0,3)左右 不能取0
	double smooth_factor = begin_smooth_factor;

	Eigen::Vector3i range(PTL.resolution + 1, PTL.resolution + 1, PTL.resolution + 1);


	//读入输入文件并输出初始模型 需要修改
	PTL.setRInner(0, -1);
	PTL.setROuter(r_outer, -1);
	PTL.setSmoothFactor(smooth_factor, -1);
	//int id = 3;
	PTL.getSkeletonFromPanetta(skeleton_path);

	if (!PTL.setSolidVolumeFractionModified(trd))
		throw std::runtime_error("Unable to initialize the requested volume fraction");
	PTL.genSolidTrussVoxelModified();			//生成体素
	PTL.genSolidTrussMeshModified();				//生成模型
	PTL.genTrussGrid();				//生成homo3d用网格
	std::string folder_name = PTL.output_path + std::to_string(id) + "-trd" + std::to_string(trd) + "-iso";
	std::filesystem::create_directories(folder_name);
	std::string obj_path = folder_name + "/" + std::to_string(id) + "-iter0.obj";
	std::string csv_path = folder_name + "/" + std::to_string(id) + "-iter0.csv";
	if (!igl::writeOBJ(obj_path, PTL.V, PTL.F))
		throw std::runtime_error("Cannot write mesh: " + obj_path);
	if (!elib::writeCSV(csv_path, PTL.voxel))
		throw std::runtime_error("Cannot write voxels: " + csv_path);

	std::cout << "初始体积分数为" << PTL.getVolumeFraction() << std::endl;


	//log文件
	std::string opt_path = folder_name + "/" + "opt.txt";
	std::ofstream file;
	file.open(opt_path);
	if (!file) throw std::runtime_error("Cannot write optimization log: " + opt_path);
	file << "-------------------- P L S M T O --------------------" << std::endl;

	//MMA初始化相关变量
	double change = 1;
	int num_iter = 0;

	//int n = 2;	//优化参数的个数 这里只有r_outer和smooth_factor
	//int m = 2;	//优化约束的个数
	int n = PTL.base_edges.size() + PTL.base_joints.size();
	int m = 3;

	double r_outer_min = 1;
	double r_outer_max = 6;
	double smooth_factor_min = 0.1;
	double smooth_factor_max = 3;

	double r_outer_movlim = 1.0;
	double smooth_factor_movlim = 0.3;

	std::vector<double> x_storage(n);
	double* x = x_storage.data();		//实际优化的所有参数组成变量
	std::vector<double> x_old_storage(n);
	double* x_old = x_old_storage.data();
	std::vector<double> f_storage(1);
	double* f = f_storage.data();		//目标量 如CH[1,1]
	std::vector<double> df_storage(n);
	double* df = df_storage.data();		//梯度
	std::vector<double> g_storage(m);
	double* g = g_storage.data();		//体积分数和各向同性约束
	std::vector<double> dg_storage(n * m);
	double* dg = dg_storage.data();	//体积分数的梯度和各向同性的梯度
	std::vector<double> x_min_storage(n);
	double* x_min = x_min_storage.data();
	std::vector<double> x_max_storage(n);
	double* x_max = x_max_storage.data();

	/*
	int x_index = 0;
	x[x_index] = PTL.r_outer;
	x_old[x_index] = PTL.r_outer;
	x_index++;

	x[x_index] = PTL.smooth_factor;
	x_old[x_index] = PTL.smooth_factor;
	x_index++;*/
	int x_index = 0;
	for (int i = 0; i < PTL.base_edges.size(); i++)
	{
		x[x_index] = PTL.base_edges[i].r_outer;
		x_old[x_index] = PTL.base_edges[i].r_outer;
		x_index++;
	}
	for (int i = 0; i < PTL.base_joints.size(); i++)
	{
		x[x_index] = PTL.base_joints[i].smooth_factor;
		x_old[x_index] = PTL.base_joints[i].smooth_factor;
		x_index++;
	}

	//MMA优化器初始化
	MMASolver mma(n, m);

	while (1)
	{
		
		homo3D::hex_element hex_e;
		setMaterial(1, 0.35, hex_e);
		init(PTL.grid, PTL.voxel, range, hex_e);
		std::vector<Eigen::VectorXd> velosity;
		auto start_all = std::chrono::high_resolution_clock::now();
		PTL.getSolidVelocitiesModified(velosity);

		
		auto start = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> duration_velo = start - start_all;
		std::cout << "开始求解lsm, velosity用时" << duration_velo.count() << "秒." << std::endl;
		solve_lsm(hex_e, velosity);	//组合出最后梯度
		auto end = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> duration = end - start;
		std::cout << "Solve_lsm运行时间: " << duration.count() << " 秒." << std::endl;

		Eigen::MatrixXd CH = hex_e.macro_constructive_tensor;

		Eigen::MatrixXd SH = CH.inverse();
		//计算dEH
		Eigen::MatrixXd dEH;
		std::vector<Eigen::MatrixXd> dSHv;
		for (int i = 0; i < hex_e.dCH.size(); i++)
		{
			Eigen::MatrixXd dshi;
			dshi.setZero(6, 6);
			dSHv.push_back(dshi);
		}
		for (int i = 0; i < n; i++) {
			dEH.setZero(6, 6);
			for (int j = 0; j < 6; j++) {
				for (int k = 0; k < 6; k++) {
					dEH(j, k) = hex_e.dCH[i](j, k);
				}
			}
			Eigen::MatrixXd dSH = -1 * SH * dEH * SH;
			for (int j = 0; j < 6; j++) {
				for (int k = 0; k < 6; k++) {
					dSHv[i](j, k) = dSH(j, k);
				}
			}
		}
		//isotropy constraint
		double ehisoii = (CH(0, 0) + CH(1, 1) + CH(2, 2)) / 3;
		double ehisoij = (CH(0, 1) + CH(0, 2) + CH(1, 2)) / 3;
		double ehisoii2 = (ehisoii - ehisoij) / 2;
		Eigen::MatrixXd EH_iso(6, 6);
		EH_iso << ehisoii, ehisoij, ehisoij, 0, 0, 0,
			ehisoij, ehisoii, ehisoij, 0, 0, 0,
			ehisoij, ehisoij, ehisoii, 0, 0, 0,
			0, 0, 0, ehisoii2, 0, 0,
			0, 0, 0, 0, ehisoii2, 0,
			0, 0, 0, 0, 0, ehisoii2;
		for (int i = 0; i < 36; i++)
		{
			Eigen::VectorXd dshi;
			dshi.setZero(n);
			hex_e.dEH.push_back(dshi);
		}
		for (int j = 0; j < 6; j++) {
			for (int k = 0; k < 6; k++) {
				for (int i = 0; i < n; i++)
				{
					hex_e.dEH[j * 6 + k](i) = hex_e.dCH[i](j, k);
				}
			}
		}
		Eigen::VectorXd dehisoii = (hex_e.dEH[0] + hex_e.dEH[7] + hex_e.dEH[14]) / 3;
		Eigen::VectorXd dehisoij = (hex_e.dEH[1] + hex_e.dEH[2] + hex_e.dEH[8]) / 3;
		//Eigen::VectorXd dehisoii2 = (dehisoii + dehisoij) / 2;		//怎么变成+了 于2024.3.31修改
		Eigen::VectorXd dehisoii2 = (dehisoii - dehisoij) / 2;
		double g_iso = 0;
		Eigen::VectorXd dgiso = Eigen::VectorXd::Zero(n);
		for (int i = 0; i < 6; i++) {
			for (int j = 0; j < 6; j++) {
				g_iso += pow(EH_iso(i, j) - CH(i, j), 2);
				if (i < 3 && j < 3 && i == j)
				{
					dgiso += 2 * (EH_iso(i, j) - CH(i, j)) / EH_iso(0, 0) * ((dehisoii - hex_e.dEH[i * 6 + j]) * EH_iso(0, 0) - dehisoii * (EH_iso(i, j) - CH(i, j))) / (EH_iso(0, 0) * EH_iso(0, 0));
				}
				else if (i < 3 && j < 3 && i != j)
				{
					dgiso += 2 * (EH_iso(i, j) - CH(i, j)) / EH_iso(0, 0) * ((dehisoij - hex_e.dEH[i * 6 + j]) * EH_iso(0, 0) - dehisoii * (EH_iso(i, j) - CH(i, j))) / (EH_iso(0, 0) * EH_iso(0, 0));
				}
				else if (i > 2 && j > 2 && i == j)
				{
					dgiso += 2 * (EH_iso(i, j) - CH(i, j)) / EH_iso(0, 0) * ((dehisoii2 - hex_e.dEH[i * 6 + j]) * EH_iso(0, 0) - dehisoii * (EH_iso(i, j) - CH(i, j))) / (EH_iso(0, 0) * EH_iso(0, 0));
				}
				else
				{
					dgiso += 2 * (-CH(i, j)) / EH_iso(0, 0) * ((-hex_e.dEH[i * 6 + j]) * EH_iso(0, 0) - dehisoii * (-CH(i, j))) / (EH_iso(0, 0) * EH_iso(0, 0));
				}
			}
		}
		g_iso /= pow(EH_iso(0, 0), 2);

		//A_U计算
		double Kv = (CH(0, 0) + CH(1, 1) + CH(2, 2)) / 9.0 + (CH(0, 1) + CH(0, 2) + CH(1, 2)) * 2.0 / 9.0;
		double Gv = (CH(0, 0) + CH(1, 1) + CH(2, 2) - CH(0, 1) - CH(0, 2) - CH(1, 2)) / 15.0 + (CH(3, 3) + CH(4, 4) + CH(5, 5)) / 5.0;
		double Kr = 1.0 / (SH(0, 0) + SH(1, 1) + SH(2, 2) + 2.0 * SH(0, 1) + 2.0 * SH(0, 2) + 2.0 * SH(1, 2));
		double Gr = 15.0 / (4.0 * SH(0, 0) + 4.0 * SH(1, 1) + 4.0 * SH(2, 2) - 4.0 * SH(0, 1) - 4.0 * SH(0, 2) - 4.0 * SH(1, 2) + 3.0 * SH(3, 3) + 3.0 * SH(4, 4) + 3.0 * SH(5, 5));
		double Au = 5.0 * Gv / Gr + Kv / Kr - 6;


		double crd = (double)PTL.voxel.sum() / (double)PTL.voxel.rows();	//当前相对密度or体积分数

		//设置目标函数
		f[0] = -1 * (CH(0, 0) + CH(1, 1) + CH(2, 2));	//在尽可能各向同性的情况下max E
		//f[0] = CH(0, 1) + CH(0, 2) + CH(1, 2) - powf(0.8, num_iter) * (CH(0, 0) + CH(1, 1) + CH(2, 2)); //minNu
		//f[0] = -1 * (CH(3, 3) + CH(4, 4) + CH(5, 5));	//maxG
		//f[0] = Au;	//iso

		//f[0] = -1 * (CH(0, 0));						//max某方向上的E
		//f[0] = -1 * (CH(3, 3));						//max某方向上的G
		//f[0] = g_iso;


		std::cout << "iter_" << num_iter << " " << "rd: " << crd << " " << "obj: " << f[0] << " " << "g_iso: " << g_iso << std::endl;
		std::cout << CH << std::endl;

		std::cout << "E: " << 1.0 / SH(0, 0) << std::endl;
		std::cout << "nu: " << -1 * SH(1, 0) / SH(0, 0) << std::endl;
		std::cout << "G: " << 1.0 / SH(3, 3) << std::endl;
		std::cout << "A^U: " << Au << std::endl;
		std::cout << "Zener ratio: " << 2 * CH(3, 3) / (CH(0, 0) - CH(0, 1)) << std::endl;
		std::cout << std::endl;

		file << "iter_" << num_iter << " " << "rd: " << crd << " " << "obj: " << f[0] << " " << "g_iso: " << g_iso << std::endl << std::endl;
		file << CH << std::endl << std::endl;

		file << "E: " << 1.0 / SH(0, 0) << std::endl;
		file << "nu: " << -1 * SH(1, 0) / SH(0, 0) << std::endl;
		file << "G: " << 1.0 / SH(3, 3) << std::endl;
		file << "A^U: " << Au << std::endl;
		file << "Zener ratio: " << 2 * CH(3, 3) / (CH(0, 0) - CH(0, 1)) << std::endl;
		file << std::endl;

		if (!((num_iter < 10 || change > 1e-2) && num_iter < max_iterations))	//终止条件,最少迭代10次 最大迭代150次?
			break;

		if (crd > trd + 0.05)	//超出体积分数约束
		{
			std::cout << "Out of target volume fraction!" << std::endl;
			file << "Out of target volume fraction!" << std::endl;
			break;
		}

		//dCH.size()=n 即优化参数个数
		for (int i = 0; i < hex_e.dCH.size(); i++) {
			df[i] = -1 * (hex_e.dCH[i](0, 0) + hex_e.dCH[i](1, 1) + hex_e.dCH[i](2, 2));
			//df[i] = hex_e.dCH[i](0, 1) + hex_e.dCH[i](0, 2) + hex_e.dCH[i](1, 2) - powf(0.8, num_iter) * (hex_e.dCH[i](0, 0) + hex_e.dCH[i](1, 1) + hex_e.dCH[i](2, 2));
			//df[i] = -1 * (hex_e.dCH[i](3, 3) + hex_e.dCH[i](4, 4) + hex_e.dCH[i](5, 5));
			//df[i] = 1.0 / 45.0 * ((CH(0, 0) + CH(1, 1) + CH(2, 2) - CH(0, 1) - CH(0, 2) - CH(1, 2)) + 3 * (CH(3, 3) + CH(4, 4) + CH(5, 5))) * (4.0 * dSHv[i](0, 0) + 4.0 * dSHv[i](1, 1) + 4.0 * dSHv[i](2, 2) - 4.0 * dSHv[i](0, 1) - 4.0 * dSHv[i](0, 2) - 4.0 * dSHv[i](1, 2) + 3.0 * dSHv[i](3, 3) + 3.0 * dSHv[i](4, 4) + 3.0 * dSHv[i](5, 5))
			//	+ 1.0 / 45.0 * ((hex_e.dCH[i](0, 0) + hex_e.dCH[i](1, 1) + hex_e.dCH[i](2, 2) - hex_e.dCH[i](0, 1) - hex_e.dCH[i](0, 2) - hex_e.dCH[i](1, 2)) + 3 * (hex_e.dCH[i](3, 3) + hex_e.dCH[i](4, 4) + hex_e.dCH[i](5, 5))) * (4.0 * SH(0, 0) + 4.0 * SH(1, 1) + 4.0 * SH(2, 2) - 4.0 * SH(0, 1) - 4.0 * SH(0, 2) - 4.0 * SH(1, 2) + 3.0 * SH(3, 3) + 3.0 * SH(4, 4) + 3.0 * SH(5, 5))
			//	+ 1.0 / 9.0 * (CH(0, 0) + CH(1, 1) + CH(2, 2) + 2.0 * (CH(0, 1) + CH(0, 2) + CH(1, 2))) * (dSHv[i](0, 0) + dSHv[i](1, 1) + dSHv[i](2, 2) + 2.0 * dSHv[i](0, 1) + 2.0 * dSHv[i](0, 2) + 2.0 * dSHv[i](1, 2))
			//	+ 1.0 / 9.0 * (hex_e.dCH[i](0, 0) + hex_e.dCH[i](1, 1) + hex_e.dCH[i](2, 2) + 2.0 * (hex_e.dCH[i](0, 1) + hex_e.dCH[i](0, 2) + hex_e.dCH[i](1, 2))) * (SH(0, 0) + SH(1, 1) + SH(2, 2) + 2.0 * SH(0, 1) + 2.0 * SH(0, 2) + 2.0 * SH(1, 2));

			//df[i] = -1 * (hex_e.dCH[i](0, 0));
			//df[i] = -1 * (hex_e.dCH[i](3, 3));
			//df[i] = dgiso[i];
			std::cout << "df[" << i << "]=" << df[i] << std::endl;

			dg[i * m] = PTL.dV[i] / trd;		//体积分数偏导
			std::cout << "dg[" << i * m << "]=" << dg[i * m] << std::endl;

			dg[i * m + 1] = -1 * PTL.dV[i] / trd;	//体积分数偏导
			std::cout << "dg[" << i * m + 1 << "]=" << dg[i * m + 1] << std::endl;

			//dg[i * m + 2] = dgiso[i];	//iso偏导 PSL版
			//std::cout << "dg[" << i * m + 2 << "]=" << dg[i * m + 2] << std::endl;

			dg[i * m + 2] = 1.0 / 45.0 * ((CH(0, 0) + CH(1, 1) + CH(2, 2) - CH(0, 1) - CH(0, 2) - CH(1, 2)) + 3 * (CH(3, 3) + CH(4, 4) + CH(5, 5))) * (4.0 * dSHv[i](0, 0) + 4.0 * dSHv[i](1, 1) + 4.0 * dSHv[i](2, 2) - 4.0 * dSHv[i](0, 1) - 4.0 * dSHv[i](0, 2) - 4.0 * dSHv[i](1, 2) + 3.0 * dSHv[i](3, 3) + 3.0 * dSHv[i](4, 4) + 3.0 * dSHv[i](5, 5))
				+ 1.0 / 45.0 * ((hex_e.dCH[i](0, 0) + hex_e.dCH[i](1, 1) + hex_e.dCH[i](2, 2) - hex_e.dCH[i](0, 1) - hex_e.dCH[i](0, 2) - hex_e.dCH[i](1, 2)) + 3 * (hex_e.dCH[i](3, 3) + hex_e.dCH[i](4, 4) + hex_e.dCH[i](5, 5))) * (4.0 * SH(0, 0) + 4.0 * SH(1, 1) + 4.0 * SH(2, 2) - 4.0 * SH(0, 1) - 4.0 * SH(0, 2) - 4.0 * SH(1, 2) + 3.0 * SH(3, 3) + 3.0 * SH(4, 4) + 3.0 * SH(5, 5))
				+ 1.0 / 9.0 * (CH(0, 0) + CH(1, 1) + CH(2, 2) + 2.0 * (CH(0, 1) + CH(0, 2) + CH(1, 2))) * (dSHv[i](0, 0) + dSHv[i](1, 1) + dSHv[i](2, 2) + 2.0 * dSHv[i](0, 1) + 2.0 * dSHv[i](0, 2) + 2.0 * dSHv[i](1, 2))
				+ 1.0 / 9.0 * (hex_e.dCH[i](0, 0) + hex_e.dCH[i](1, 1) + hex_e.dCH[i](2, 2) + 2.0 * (hex_e.dCH[i](0, 1) + hex_e.dCH[i](0, 2) + hex_e.dCH[i](1, 2))) * (SH(0, 0) + SH(1, 1) + SH(2, 2) + 2.0 * SH(0, 1) + 2.0 * SH(0, 2) + 2.0 * SH(1, 2));

		}
		g[0] = crd / trd - 1;
		g[1] = 0.98 - crd / trd;
		//g[2] = g_iso - 1e-5;
		g[2] = Au - 1e-5;
		//std::cout << "g[0]=" << g[0] << ",g[1]=" << g[1] << std::endl;
		std::cout << "g[0]=" << g[0] << ",g[1]=" << g[1] << ",g[2]=" << g[2] << std::endl;

		//设置每次步进的限制范围
		x_index = 0;
		for (int i = 0; i < PTL.base_edges.size(); i++)
		{
			x_max[x_index] = std::min(r_outer_max, x[x_index] + r_outer_movlim);
			x_min[x_index] = std::max(r_outer_min, x[x_index] - r_outer_movlim);
			x_index++;
		}
		for (int i = 0; i < PTL.base_joints.size(); i++)
		{
			x_max[x_index] = std::min(smooth_factor_max, x[x_index] + smooth_factor_movlim);
			x_min[x_index] = std::max(smooth_factor_min, x[x_index] - smooth_factor_movlim);
			x_index++;
		}


		mma.Update(x, df, g, dg, x_min, x_max);		//求解



		//更新参数
		x_index = 0;
		for (int i = 0; i < PTL.base_edges.size(); i++)
		{
			std::cout << "R_outer[" << i << "]:" << x_old[x_index] << " -> " << x[x_index] << std::endl;
			file << "R_outer[" << i << "]:" << x_old[x_index] << " -> " << x[x_index] << std::endl;
			PTL.setROuter(x[x_index], i);
			x_index++;
		}
		for (int i = 0; i < PTL.base_joints.size(); i++)
		{
			std::cout << "smooth_factor[" << i << "]:" << x_old[x_index] << " -> " << x[x_index] << std::endl;
			file << "smooth_factor[" << i << "]:" << x_old[x_index] << " -> " << x[x_index] << std::endl;
			PTL.setSmoothFactor(x[x_index], i);
			x_index++;
		}

		double maxdiff = 1e-10;
		for (int i = 0; i < n; i++) {
			double diff = std::abs(x[i] - x_old[i]);
			maxdiff = std::max(maxdiff, diff);
			x_old[i] = x[i];
		}
		change = maxdiff;

		std::cout << "change:" << change << std::endl;
		std::cout << std::endl;
		file << "change:" << change << std::endl;
		file << std::endl;

		//重新生成
		num_iter++;
		PTL.genSolidTrussVoxelModified();			//生成体素
		PTL.genSolidTrussMeshModified();				//生成模型
		PTL.genTrussGrid();				//生成homo3d用网格
		obj_path = folder_name + "/" + std::to_string(id) + "-iter" + std::to_string(num_iter) + ".obj";
		csv_path = folder_name + "/" + std::to_string(id) + "-iter" + std::to_string(num_iter) + ".csv";
		if (!igl::writeOBJ(obj_path, PTL.V, PTL.F))
			throw std::runtime_error("Cannot write mesh: " + obj_path);
		if (!elib::writeCSV(csv_path, PTL.voxel))
			throw std::runtime_error("Cannot write voxels: " + csv_path);
		auto end_all = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> duration_all = end_all - start_all;
		std::cout << "迭代一次运行时间: " << duration_all.count() << " 秒" << std::endl;
	}

	file.close();
}

void plsmtoSolidModifiedETH(Truss& PTL, double trd, int id, double begin_smooth_factor, const std::string& skeleton_path, int max_iterations)
{
	//默认参数
	double unit_size = PTL.cell_length;	//一个单胞的尺寸

	double relative_r_outer = 0.04;			//相对半径,然后再*unit_size
	double relative_smooth_factor = 0.01;

	double r_outer = unit_size * relative_r_outer; //范围1.5-4左右
	//double smooth_factor = unit_size * relative_smooth_factor;//(0,3)左右 不能取0
	double smooth_factor = begin_smooth_factor;

	Eigen::Vector3i range(PTL.resolution + 1, PTL.resolution + 1, PTL.resolution + 1);


	//读入输入文件并输出初始模型 需要修改
	PTL.setRInner(0, -1);
	PTL.setROuter(r_outer, -1);
	PTL.setSmoothFactor(smooth_factor, -1);
	//int id = 3;
	PTL.getSkeletonFromETH(skeleton_path);

	if (!PTL.setSolidVolumeFractionModifiedETH(trd))
		throw std::runtime_error("Unable to initialize the requested volume fraction");
	PTL.genSolidTrussVoxelModifiedETH();			//生成体素
	PTL.genSolidTrussMeshModifiedETH();				//生成模型
	PTL.genTrussGrid();				//生成homo3d用网格
	std::string folder_name = PTL.output_path + std::to_string(id) + "-ETHtrd" + std::to_string(trd) + "-iso";
	std::filesystem::create_directories(folder_name);
	std::string obj_path = folder_name + "/" + std::to_string(id) + "-iter0.obj";
	std::string csv_path = folder_name + "/" + std::to_string(id) + "-iter0.csv";
	if (!igl::writeOBJ(obj_path, PTL.V, PTL.F))
		throw std::runtime_error("Cannot write mesh: " + obj_path);
	if (!elib::writeCSV(csv_path, PTL.voxel))
		throw std::runtime_error("Cannot write voxels: " + csv_path);

	std::cout << "初始体积分数为" << PTL.getVolumeFraction() << std::endl;


	//log文件
	std::string opt_path = folder_name + "/" + "opt.txt";
	std::ofstream file;
	file.open(opt_path);
	if (!file) throw std::runtime_error("Cannot write optimization log: " + opt_path);
	file << "-------------------- P L S M T O --------------------" << std::endl;

	//MMA初始化相关变量
	double change = 1;
	int num_iter = 0;

	//int n = 2;	//优化参数的个数 这里只有r_outer和smooth_factor
	//int m = 2;	//优化约束的个数
	int n = PTL.base_edges.size() + PTL.base_joints.size();
	int m = 2;

	double r_outer_min = 1.2;
	double r_outer_max = 6;
	double smooth_factor_min = 0.1;
	double smooth_factor_max = 0.5;

	double r_outer_movlim = 1.0;
	double smooth_factor_movlim = 0.3;

	std::vector<double> x_storage(n);
	double* x = x_storage.data();		//实际优化的所有参数组成变量
	std::vector<double> x_old_storage(n);
	double* x_old = x_old_storage.data();
	std::vector<double> f_storage(1);
	double* f = f_storage.data();		//目标量 如CH[1,1]
	std::vector<double> df_storage(n);
	double* df = df_storage.data();		//梯度
	std::vector<double> g_storage(m);
	double* g = g_storage.data();		//体积分数和各向同性约束
	std::vector<double> dg_storage(n * m);
	double* dg = dg_storage.data();	//体积分数的梯度和各向同性的梯度
	std::vector<double> x_min_storage(n);
	double* x_min = x_min_storage.data();
	std::vector<double> x_max_storage(n);
	double* x_max = x_max_storage.data();

	/*
	int x_index = 0;
	x[x_index] = PTL.r_outer;
	x_old[x_index] = PTL.r_outer;
	x_index++;

	x[x_index] = PTL.smooth_factor;
	x_old[x_index] = PTL.smooth_factor;
	x_index++;*/
	int x_index = 0;
	for (int i = 0; i < PTL.base_edges.size(); i++)
	{
		x[x_index] = PTL.base_edges[i].r_outer;
		x_old[x_index] = PTL.base_edges[i].r_outer;
		x_index++;
	}
	for (int i = 0; i < PTL.base_joints.size(); i++)
	{
		x[x_index] = PTL.base_joints[i].smooth_factor;
		x_old[x_index] = PTL.base_joints[i].smooth_factor;
		x_index++;
	}

	//MMA优化器初始化
	MMASolver mma(n, m);

	while (1)
	{
		homo3D::hex_element hex_e;
		setMaterial(1, 0.35, hex_e);
		init(PTL.grid, PTL.voxel, range, hex_e);
		std::vector<Eigen::VectorXd> velosity;

		PTL.getSolidVelocitiesModifiedETH(velosity);

		solve_lsm(hex_e, velosity);	//组合出最后梯度

		Eigen::MatrixXd CH = hex_e.macro_constructive_tensor;

		Eigen::MatrixXd SH = CH.inverse();
		//计算dEH
		Eigen::MatrixXd dEH;
		std::vector<Eigen::MatrixXd> dSHv;
		for (int i = 0; i < hex_e.dCH.size(); i++)
		{
			Eigen::MatrixXd dshi;
			dshi.setZero(6, 6);
			dSHv.push_back(dshi);
		}
		for (int i = 0; i < n; i++) {
			dEH.setZero(6, 6);
			for (int j = 0; j < 6; j++) {
				for (int k = 0; k < 6; k++) {
					dEH(j, k) = hex_e.dCH[i](j, k);
				}
			}
			Eigen::MatrixXd dSH = -1 * SH * dEH * SH;
			for (int j = 0; j < 6; j++) {
				for (int k = 0; k < 6; k++) {
					dSHv[i](j, k) = dSH(j, k);
				}
			}
		}
		//isotropy constraint
		double ehisoii = (CH(0, 0) + CH(1, 1) + CH(2, 2)) / 3;
		double ehisoij = (CH(0, 1) + CH(0, 2) + CH(1, 2)) / 3;
		double ehisoii2 = (ehisoii - ehisoij) / 2;
		Eigen::MatrixXd EH_iso(6, 6);
		EH_iso << ehisoii, ehisoij, ehisoij, 0, 0, 0,
			ehisoij, ehisoii, ehisoij, 0, 0, 0,
			ehisoij, ehisoij, ehisoii, 0, 0, 0,
			0, 0, 0, ehisoii2, 0, 0,
			0, 0, 0, 0, ehisoii2, 0,
			0, 0, 0, 0, 0, ehisoii2;
		for (int i = 0; i < 36; i++)
		{
			Eigen::VectorXd dshi;
			dshi.setZero(n);
			hex_e.dEH.push_back(dshi);
		}
		for (int j = 0; j < 6; j++) {
			for (int k = 0; k < 6; k++) {
				for (int i = 0; i < n; i++)
				{
					hex_e.dEH[j * 6 + k](i) = hex_e.dCH[i](j, k);
				}
			}
		}
		Eigen::VectorXd dehisoii = (hex_e.dEH[0] + hex_e.dEH[7] + hex_e.dEH[14]) / 3;
		Eigen::VectorXd dehisoij = (hex_e.dEH[1] + hex_e.dEH[2] + hex_e.dEH[8]) / 3;
		//Eigen::VectorXd dehisoii2 = (dehisoii + dehisoij) / 2;		//怎么变成+了 于2024.3.31修改
		Eigen::VectorXd dehisoii2 = (dehisoii - dehisoij) / 2;
		double g_iso = 0;
		Eigen::VectorXd dgiso = Eigen::VectorXd::Zero(n);
		for (int i = 0; i < 6; i++) {
			for (int j = 0; j < 6; j++) {
				g_iso += pow(EH_iso(i, j) - CH(i, j), 2);
				if (i < 3 && j < 3 && i == j)
				{
					dgiso += 2 * (EH_iso(i, j) - CH(i, j)) / EH_iso(0, 0) * ((dehisoii - hex_e.dEH[i * 6 + j]) * EH_iso(0, 0) - dehisoii * (EH_iso(i, j) - CH(i, j))) / (EH_iso(0, 0) * EH_iso(0, 0));
				}
				else if (i < 3 && j < 3 && i != j)
				{
					dgiso += 2 * (EH_iso(i, j) - CH(i, j)) / EH_iso(0, 0) * ((dehisoij - hex_e.dEH[i * 6 + j]) * EH_iso(0, 0) - dehisoii * (EH_iso(i, j) - CH(i, j))) / (EH_iso(0, 0) * EH_iso(0, 0));
				}
				else if (i > 2 && j > 2 && i == j)
				{
					dgiso += 2 * (EH_iso(i, j) - CH(i, j)) / EH_iso(0, 0) * ((dehisoii2 - hex_e.dEH[i * 6 + j]) * EH_iso(0, 0) - dehisoii * (EH_iso(i, j) - CH(i, j))) / (EH_iso(0, 0) * EH_iso(0, 0));
				}
				else
				{
					dgiso += 2 * (-CH(i, j)) / EH_iso(0, 0) * ((-hex_e.dEH[i * 6 + j]) * EH_iso(0, 0) - dehisoii * (-CH(i, j))) / (EH_iso(0, 0) * EH_iso(0, 0));
				}
			}
		}
		g_iso /= pow(EH_iso(0, 0), 2);

		//A_U计算
		double Kv = (CH(0, 0) + CH(1, 1) + CH(2, 2)) / 9.0 + (CH(0, 1) + CH(0, 2) + CH(1, 2)) * 2.0 / 9.0;
		double Gv = (CH(0, 0) + CH(1, 1) + CH(2, 2) - CH(0, 1) - CH(0, 2) - CH(1, 2)) / 15.0 + (CH(3, 3) + CH(4, 4) + CH(5, 5)) / 5.0;
		double Kr = 1.0 / (SH(0, 0) + SH(1, 1) + SH(2, 2) + 2.0 * SH(0, 1) + 2.0 * SH(0, 2) + 2.0 * SH(1, 2));
		double Gr = 15.0 / (4.0 * SH(0, 0) + 4.0 * SH(1, 1) + 4.0 * SH(2, 2) - 4.0 * SH(0, 1) - 4.0 * SH(0, 2) - 4.0 * SH(1, 2) + 3.0 * SH(3, 3) + 3.0 * SH(4, 4) + 3.0 * SH(5, 5));
		double Au = 5.0 * Gv / Gr + Kv / Kr - 6;


		double crd = (double)PTL.voxel.sum() / (double)PTL.voxel.rows();	//当前相对密度or体积分数

		//设置目标函数
		//f[0] = -1 * (CH(0, 0) + CH(1, 1) + CH(2, 2));	//在尽可能各向同性的情况下max E
		//f[0] = CH(0, 1) + CH(0, 2) + CH(1, 2) - powf(0.5, num_iter) * (CH(0, 0) + CH(1, 1) + CH(2, 2)); //minNu
		//f[0] = -1 * (CH(3, 3) + CH(4, 4) + CH(5, 5));	//maxG
		f[0] = Au;	//iso

		//f[0] = -1 * (CH(0, 0));						//max某方向上的E
		//f[0] = -1 * (CH(3, 3));						//max某方向上的G
		//f[0] = g_iso;


		std::cout << "iter_" << num_iter << " " << "rd: " << crd << " " << "obj: " << f[0] << " " << "g_iso: " << g_iso << std::endl;
		std::cout << CH << std::endl;

		std::cout << "E: " << 1.0 / SH(0, 0) << std::endl;
		std::cout << "nu: " << -1 * SH(1, 0) / SH(0, 0) << std::endl;
		std::cout << "G: " << 1.0 / SH(3, 3) << std::endl;
		std::cout << "A^U: " << Au << std::endl;
		std::cout << "Zener ratio: " << 2 * CH(3, 3) / (CH(0, 0) - CH(0, 1)) << std::endl;
		std::cout << std::endl;

		file << "iter_" << num_iter << " " << "rd: " << crd << " " << "obj: " << f[0] << " " << "g_iso: " << g_iso << std::endl << std::endl;
		file << CH << std::endl << std::endl;

		file << "E: " << 1.0 / SH(0, 0) << std::endl;
		file << "nu: " << -1 * SH(1, 0) / SH(0, 0) << std::endl;
		file << "G: " << 1.0 / SH(3, 3) << std::endl;
		file << "A^U: " << Au << std::endl;
		file << "Zener ratio: " << 2 * CH(3, 3) / (CH(0, 0) - CH(0, 1)) << std::endl;
		file << std::endl;

		if (!((num_iter < 10 || change > 1e-4) && num_iter < max_iterations))	//终止条件,最少迭代10次 最大迭代150次?
			break;

		if (crd > trd + 0.05)	//超出体积分数约束
		{
			std::cout << "Out of target volume fraction!" << std::endl;
			file << "Out of target volume fraction!" << std::endl;
			break;
		}

		//dCH.size()=n 即优化参数个数
		for (int i = 0; i < hex_e.dCH.size(); i++) {
			//df[i] = -1 * (hex_e.dCH[i](0, 0) + hex_e.dCH[i](1, 1) + hex_e.dCH[i](2, 2));
			//df[i] = hex_e.dCH[i](0, 1) + hex_e.dCH[i](0, 2) + hex_e.dCH[i](1, 2) - powf(0.5, num_iter) * (hex_e.dCH[i](0, 0) + hex_e.dCH[i](1, 1) + hex_e.dCH[i](2, 2));
			//df[i] = -1 * (hex_e.dCH[i](3, 3) + hex_e.dCH[i](4, 4) + hex_e.dCH[i](5, 5));
			df[i] = 1.0 / 45.0 * ((CH(0, 0) + CH(1, 1) + CH(2, 2) - CH(0, 1) - CH(0, 2) - CH(1, 2)) + 3 * (CH(3, 3) + CH(4, 4) + CH(5, 5))) * (4.0 * dSHv[i](0, 0) + 4.0 * dSHv[i](1, 1) + 4.0 * dSHv[i](2, 2) - 4.0 * dSHv[i](0, 1) - 4.0 * dSHv[i](0, 2) - 4.0 * dSHv[i](1, 2) + 3.0 * dSHv[i](3, 3) + 3.0 * dSHv[i](4, 4) + 3.0 * dSHv[i](5, 5))
				+ 1.0 / 45.0 * ((hex_e.dCH[i](0, 0) + hex_e.dCH[i](1, 1) + hex_e.dCH[i](2, 2) - hex_e.dCH[i](0, 1) - hex_e.dCH[i](0, 2) - hex_e.dCH[i](1, 2)) + 3 * (hex_e.dCH[i](3, 3) + hex_e.dCH[i](4, 4) + hex_e.dCH[i](5, 5))) * (4.0 * SH(0, 0) + 4.0 * SH(1, 1) + 4.0 * SH(2, 2) - 4.0 * SH(0, 1) - 4.0 * SH(0, 2) - 4.0 * SH(1, 2) + 3.0 * SH(3, 3) + 3.0 * SH(4, 4) + 3.0 * SH(5, 5))
				+ 1.0 / 9.0 * (CH(0, 0) + CH(1, 1) + CH(2, 2) + 2.0 * (CH(0, 1) + CH(0, 2) + CH(1, 2))) * (dSHv[i](0, 0) + dSHv[i](1, 1) + dSHv[i](2, 2) + 2.0 * dSHv[i](0, 1) + 2.0 * dSHv[i](0, 2) + 2.0 * dSHv[i](1, 2))
				+ 1.0 / 9.0 * (hex_e.dCH[i](0, 0) + hex_e.dCH[i](1, 1) + hex_e.dCH[i](2, 2) + 2.0 * (hex_e.dCH[i](0, 1) + hex_e.dCH[i](0, 2) + hex_e.dCH[i](1, 2))) * (SH(0, 0) + SH(1, 1) + SH(2, 2) + 2.0 * SH(0, 1) + 2.0 * SH(0, 2) + 2.0 * SH(1, 2));
			
			//df[i] = -1 * (hex_e.dCH[i](0, 0));
			//df[i] = -1 * (hex_e.dCH[i](3, 3));
			//df[i] = dgiso[i];
			std::cout << "df[" << i << "]=" << df[i] << std::endl;

			dg[i * m] = PTL.dV[i] / trd;		//体积分数偏导
			std::cout << "dg[" << i * m << "]=" << dg[i * m] << std::endl;

			dg[i * m + 1] = -1 * PTL.dV[i] / trd;	//体积分数偏导
			std::cout << "dg[" << i * m + 1 << "]=" << dg[i * m + 1] << std::endl;

			//dg[i * m + 2] = dgiso[i];	//iso偏导 PSL版
			//std::cout << "dg[" << i * m + 2 << "]=" << dg[i * m + 2] << std::endl;

			/*dg[i * m + 2] = 1.0 / 45.0 * ((CH(0, 0) + CH(1, 1) + CH(2, 2) - CH(0, 1) - CH(0, 2) - CH(1, 2)) + 3 * (CH(3, 3) + CH(4, 4) + CH(5, 5))) * (4.0 * dSHv[i](0, 0) + 4.0 * dSHv[i](1, 1) + 4.0 * dSHv[i](2, 2) - 4.0 * dSHv[i](0, 1) - 4.0 * dSHv[i](0, 2) - 4.0 * dSHv[i](1, 2) + 3.0 * dSHv[i](3, 3) + 3.0 * dSHv[i](4, 4) + 3.0 * dSHv[i](5, 5))
				+ 1.0 / 45.0 * ((hex_e.dCH[i](0, 0) + hex_e.dCH[i](1, 1) + hex_e.dCH[i](2, 2) - hex_e.dCH[i](0, 1) - hex_e.dCH[i](0, 2) - hex_e.dCH[i](1, 2)) + 3 * (hex_e.dCH[i](3, 3) + hex_e.dCH[i](4, 4) + hex_e.dCH[i](5, 5))) * (4.0 * SH(0, 0) + 4.0 * SH(1, 1) + 4.0 * SH(2, 2) - 4.0 * SH(0, 1) - 4.0 * SH(0, 2) - 4.0 * SH(1, 2) + 3.0 * SH(3, 3) + 3.0 * SH(4, 4) + 3.0 * SH(5, 5))
				+ 1.0 / 9.0 * (CH(0, 0) + CH(1, 1) + CH(2, 2) + 2.0 * (CH(0, 1) + CH(0, 2) + CH(1, 2))) * (dSHv[i](0, 0) + dSHv[i](1, 1) + dSHv[i](2, 2) + 2.0 * dSHv[i](0, 1) + 2.0 * dSHv[i](0, 2) + 2.0 * dSHv[i](1, 2))
				+ 1.0 / 9.0 * (hex_e.dCH[i](0, 0) + hex_e.dCH[i](1, 1) + hex_e.dCH[i](2, 2) + 2.0 * (hex_e.dCH[i](0, 1) + hex_e.dCH[i](0, 2) + hex_e.dCH[i](1, 2))) * (SH(0, 0) + SH(1, 1) + SH(2, 2) + 2.0 * SH(0, 1) + 2.0 * SH(0, 2) + 2.0 * SH(1, 2));
			*/
		}
		g[0] = crd / trd - 1;
		g[1] = 0.98 - crd / trd;
		//g[2] = g_iso - 1e-5;
		//g[2] = Au - 1e-5;
		std::cout << "g[0]=" << g[0] << ",g[1]=" << g[1] << std::endl;
		//std::cout << "g[0]=" << g[0] << ",g[1]=" << g[1] << ",g[2]=" << g[2] << std::endl;

		//设置每次步进的限制范围
		x_index = 0;
		for (int i = 0; i < PTL.base_edges.size(); i++)
		{
			x_max[x_index] = std::min(r_outer_max, x[x_index] + r_outer_movlim);
			x_min[x_index] = std::max(r_outer_min, x[x_index] - r_outer_movlim);
			x_index++;
		}
		for (int i = 0; i < PTL.base_joints.size(); i++)
		{
			x_max[x_index] = std::min(smooth_factor_max, x[x_index] + smooth_factor_movlim);
			x_min[x_index] = std::max(smooth_factor_min, x[x_index] - smooth_factor_movlim);
			x_index++;
		}


		mma.Update(x, df, g, dg, x_min, x_max);		//求解



		//更新参数
		x_index = 0;
		for (int i = 0; i < PTL.base_edges.size(); i++)
		{
			std::cout << "R_outer[" << i << "]:" << x_old[x_index] << " -> " << x[x_index] << std::endl;
			file << "R_outer[" << i << "]:" << x_old[x_index] << " -> " << x[x_index] << std::endl;
			PTL.setROuter(x[x_index], i);
			x_index++;
		}
		for (int i = 0; i < PTL.base_joints.size(); i++)
		{
			std::cout << "smooth_factor[" << i << "]:" << x_old[x_index] << " -> " << x[x_index] << std::endl;
			file << "smooth_factor[" << i << "]:" << x_old[x_index] << " -> " << x[x_index] << std::endl;
			PTL.setSmoothFactor(x[x_index], i);
			x_index++;
		}

		double maxdiff = 1e-10;
		for (int i = 0; i < n; i++) {
			double diff = std::abs(x[i] - x_old[i]);
			maxdiff = std::max(maxdiff, diff);
			x_old[i] = x[i];
		}
		change = maxdiff;

		std::cout << "change:" << change << std::endl;
		std::cout << std::endl;
		file << "change:" << change << std::endl;
		file << std::endl;

		//重新生成
		num_iter++;
		PTL.genSolidTrussVoxelModifiedETH();			//生成体素
		PTL.genSolidTrussMeshModifiedETH();				//生成模型
		PTL.genTrussGrid();				//生成homo3d用网格
		obj_path = folder_name + "/" + std::to_string(id) + "-iter" + std::to_string(num_iter) + ".obj";
		csv_path = folder_name + "/" + std::to_string(id) + "-iter" + std::to_string(num_iter) + ".csv";
		if (!igl::writeOBJ(obj_path, PTL.V, PTL.F))
			throw std::runtime_error("Cannot write mesh: " + obj_path);
		if (!elib::writeCSV(csv_path, PTL.voxel))
			throw std::runtime_error("Cannot write voxels: " + csv_path);
	}

	file.close();
}
