#pragma once
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <vector>
#include <iostream>
#include <sstream>
#include <string>
#include <chrono>
#include "write_csv.h"
#include "read_csv.h"

const long double petl_pi = 3.14159265358979323L;

struct Edge
{
	Eigen::MatrixXd v;				//2*3矩阵
	double r_outer;					//外半径
	double r_inner;					//内半径,暂时用不上
	int r_index;						//记录这条边的半径与哪条base_edge相关联
};


struct Joint
{
	Eigen::Vector3d center;				//中心
	double smooth_factor;				//平滑因子
	std::vector<int> associated_edges;	//相邻边,值为edges的坐标
	int smf_index;						//记录这个交点的smf与哪个base_joint相关联
	bool is_same_side;					//记录该交点的所有相邻边是否在空间同一侧
};

class Truss
{
public:
	//基础数据
	int num_edges = 0;
	int num_nodes = 0;
	//int num_control_points = 0;	
	std::vector<Edge> base_edges;			//子单元的骨架
	std::vector<Joint> base_joints;			//子单元的关节
	std::vector<Edge> edges;					//所有边
	std::vector<Joint> nodes;				//所有顶点
	//std::vector<std::vector<Eigen::MatrixXd>> topo;	//拓扑结构
	//Eigen::MatrixXd control_points;	//控制点
	double cell_length = 50;			//模型包围盒大小
	double eps = 1e-5;				//极小值


	//控制参数
	double r_inner = 0;			//杆的内半径,用于初始化
	double r_outer = 3;			//杆的外半径,用于初始化
	double smooth_factor = 3;	//关节处平滑因子

	//double rotation_alpha = 0;	//旋转参数alpha
	//double rotation_beta = 0;	//旋转参数beta
	//double slide_dis = 0;		//滑动距离

	//生成mesh需要使用的数据
	Eigen::MatrixXd V;			//顶点
	Eigen::MatrixXi F;			//面片
	int resolution = 64;			//体素化精度
	int resolution_mesh = 64;	//Marching Cubes重建精度
	Eigen::MatrixXd grid;		//存放体素坐标
	Eigen::VectorXi voxel;		//存放体素的值
	std::string input_path;		//输入骨架的文件夹路径
	std::string output_path;		//输出模型的文件夹路径

	//优化相关
	std::vector<double> dV;		//体积分数的变化值
	std::vector<double> dH;		//好像用不上
	double tvol = 0;				//也是体积分数	
	//Eigen::VectorXi vol;		//体素的值,好像跟voxel是一样的
	Eigen::VectorXd outer_sdf;	//外部距离场
	//Eigen::VectorXd inner_sdf;	//内部距离场,镂空杆需要分开存放来计算距离场的差

	//函数
	Truss();
	~Truss();
	void setROuter(double r, int index)	//给base_edges[index]相关联的杆赋值R_outer
	{
		if (index == -1)		//所有边全设置成r
		{
			r_outer = r;
			for (int i = 0; i < base_edges.size(); i++)
				base_edges[i].r_outer = r_outer;
			for (int i = 0; i < edges.size(); i++)
				edges[i].r_outer = r_outer;
			return;
		}
		base_edges[index].r_outer = r;
		for (int i = 0; i < edges.size(); i++)
		{
			if (edges[i].r_index == index)
				edges[i].r_outer = r;
		}

	}
	void setRInner(double r, int index)		//给base_edges[index]相关联的杆赋值R_inner
	{
		if (index == -1)		//所有边全设置成r
		{
			r_inner = r;
			for (int i = 0; i < base_edges.size(); i++)
				base_edges[i].r_inner = r_inner;
			for (int i = 0; i < edges.size(); i++)
				edges[i].r_inner = r_inner;
			return;
		}
		base_edges[index].r_inner = r;
		for (int i = 0; i < edges.size(); i++)
		{
			if (edges[i].r_index == index)
				edges[i].r_inner = r;
		}

	}
	void setSmoothFactor(double s, int index)	//给base_joints[index]相关联的杆赋值smf
	{
		if (index == -1) 	//所有顶点全设置s
		{
			smooth_factor = s;
			for (int i = 0; i < base_joints.size(); i++)
				base_joints[i].smooth_factor = s;
			for (int i = 0; i < nodes.size(); i++)
				nodes[i].smooth_factor = s;
			return;
		}
		base_joints[index].smooth_factor = s;
		for (int i = 0; i < nodes.size(); i++)
		{
			if (nodes[i].smf_index == index)
				nodes[i].smooth_factor = s;
		}
	}
	void setCellLength(double c)	//设置包围盒大小
	{
		cell_length = c;
	}

	double getVolumeFraction()
	{
		return (double)voxel.sum() / (double)voxel.rows();
	}
	double getROuter(int index)
	{
		//检查index是否越界
		if (base_edges.size() <= index)
		{
			return -1;
		}
		return base_edges[index].r_outer;
	}
	double getRInner(int index)
	{
		//检查index是否越界
		if (base_edges.size() <= index)
		{
			return -1;
		}
		return base_edges[index].r_inner;
	}
	double getSmoothFactor(int index)
	{
		if (base_joints.size() <= index)
		{
			return -1;
		}
		return base_joints[index].smooth_factor;
	}


	bool loadConfig(std::string LoadPath);	//读取输入信息config

	//bool isEquation(Eigen::Vector3d line1, Eigen::Vector3d line2);
	//bool inBaseTetrahedral(Eigen::Vector3d v);			//该点是否在初始1/48四面体中,在.cpp里声明定义
	//void getPanettaSymmetryPoint(Eigen::Vector3d &v);	//通过1/48对称规则得到初始四面体中的对应点,在.cpp里声明定义
	//bool edgesAreSameSide(Eigen::Vector3d joint_center,Eigen::MatrixXd points);	//检查joint的相邻边是否在空间同一侧
	bool edgesCanMerge(Edge e1, Edge e2, Edge& e, Eigen::Vector3d& e_joint);		//两条边可以合并
	void getTopology(std::vector<Edge> edge_list);
	void getSkeletonFromPanetta(std::string path);	//读取Panetta的骨架

	//bool inBaseCube(Eigen::Vector3d v);	//该点是否在初始1/8立方体中,在.cpp里声明定义
	//void getETHSymmetryPoint(Eigen::Vector3d &v);	//通过1/8对称规则得到初始四面体中的对应点,在.cpp里声明定义
	void getTopologyETH(std::vector<Edge> edge_list);
	void getSkeletonFromETH(std::string path);	//读取ETH的骨架

	//实心杆-1/48规则
	bool setSolidVolumeFraction(double vof);		//通过调整半径r_outer来设置体积分数
	void genSolidTrussVoxel();					//生成杆的体素文件,得到voxel
	void genSolidTrussMesh();					//生成杆的模型,得到顶点V和面片F,与上面不同的是构造距离场时需要比模型边界多一圈
	void genTilingTrussMesh(int size_x, int size_y, int size_z);		//平铺size_x,size_y,size_z

	//实心杆,去膨胀
	bool setSolidVolumeFractionModified(double vof);		//通过调整半径r_outer来设置体积分数
	void genSolidTrussVoxelModified();
	void genSolidTrussMeshModified();
	void genTilingTrussMeshModified(int size_x, int size_y, int size_z);		//平铺size_x,size_y,size_z

	//融合,对所有edge应用统一smf-1/48规则
	bool setFusionVolumeFraction(double vof);	//通过调整半径r_outer来设置体积分数
	void genFusionTrussVoxel();					//生成杆的体素文件,得到voxel
	void genFusionTrussMesh();					//生成杆的模型,得到顶点V和面片F,与上面不同的是构造距离场时需要比模型边界多一圈

	//实心杆-1/8规则
	bool setSolidVolumeFractionETH(double vof);		//通过调整半径r_outer来设置体积分数
	void genSolidTrussVoxelETH();					//生成杆的体素文件,得到voxel
	void genSolidTrussMeshETH();					//生成杆的模型,得到顶点V和面片F,与上面不同的是构造距离场时需要比模型边界多一圈

	//实心杆,去膨胀
	bool setSolidVolumeFractionModifiedETH(double vof);		//通过调整半径r_outer来设置体积分数
	void genSolidTrussVoxelModifiedETH();					//生成杆的体素文件,得到voxel
	void genSolidTrussMeshModifiedETH();					//生成杆的模型,得到顶点V和面片F,与上面不同的是构造距离场时需要比模型边界多一圈
	void genTilingTrussMeshModifiedETH(int size_x, int size_y, int size_z);		//平铺size_x,size_y,size_z

	//普通杆,不使用smf,直接算到骨架的最短距离-1/8规则
	bool setNormalVolumeFractionETH(double vof);	//通过调整半径r_outer来设置体积分数
	void genNormalTrussVoxelETH();					//生成杆的体素文件,得到voxel
	void genNormalTrussMeshETH();					//生成杆的模型,得到顶点V和面片F,与上面不同的是构造距离场时需要比模型边界多一圈

	//镂空杆
	//void setHollowVolumeFraction(double vof);	//固定r_inner,通过调整半径r_outer来设置体积分数
	//void genHollowTrussVoxel();					//生成杆的体素文件,得到voxel
	//void genHollowTrussMesh();					//生成杆的模型,得到顶点V和面片F,与上面不同的是构造距离场时需要比模型边界多一圈

	//优化相关
	bool isInnerVoxel(int index);				//体素index是否是内部体素,通过邻居体素是否值为1来判断
	void getTargetCoors(Eigen::MatrixXd& tgrid, Eigen::VectorXi& mapping);	//生成值全为1的体素坐标,mapping[全部体素坐标i]=[值为1体素坐标j]
	void genTrussGrid();													//生成grid网格为分辨率+1

	void genSolidTargetSdf(Eigen::MatrixXd& tgrid);								//生成目标距离场
	void getSolidVelocities(std::vector<Eigen::VectorXd>& velosities);			//获取速度

	void genSolidTargetSdfModified(Eigen::MatrixXd& tgrid);								//生成目标距离场
	void getSolidVelocitiesModified(std::vector<Eigen::VectorXd>& velosities);			//获取速度


	void genSolidTargetSdfETH(Eigen::MatrixXd& tgrid);								//生成目标距离场
	void getSolidVelocitiesETH(std::vector<Eigen::VectorXd>& velosities);			//获取速度

	void genSolidTargetSdfModifiedETH(Eigen::MatrixXd& tgrid);								//生成目标距离场
	void getSolidVelocitiesModifiedETH(std::vector<Eigen::VectorXd>& velosities);			//获取速度

	//void genHollowTargetSdf(Eigen::MatrixXd& tgrid);								//生成目标距离场
	//void getHollowVelocities(std::vector<Eigen::VectorXd>& velosities);			//获取速度

	
	
};
