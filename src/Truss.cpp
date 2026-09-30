#include "petl/Truss.h"
#include <cassert>
#include <unsupported/Eigen/NNLS>
#include <igl/copyleft/marching_cubes.h>
#include <stdexcept>

Truss::Truss()
{

}

Truss::~Truss()
{

}

//检查joint的相邻边是否在空间同一侧
//joint_center:1*3向量
//points:n*3矩阵
bool edgesAreSameSide(Eigen::Vector3d joint_center, Eigen::MatrixXd points)
{
	if (points.rows() <= 3)
		return true;
	//计算边向量
	Eigen::MatrixXd vectors = points.rowwise() - joint_center.transpose();
	//考虑所有向量对组合
	int numVectors = vectors.rows();
	for (int i = 0; i < numVectors - 1; ++i) {
		for (int j = i + 1; j < numVectors; ++j) {
			// 计算法向量
			Eigen::Vector3d vi = vectors.row(i);
			Eigen::Vector3d vj = vectors.row(j);
			Eigen::Vector3d n = vi.cross(vj);

			// 如果法向量接近零向量，则ei和ej共线，跳过
			if (n.norm() < 1e-10) {
				continue;
			}

			// 计算所有向量与法向量的点积
			Eigen::VectorXd signs = vectors * n;

			// 如果所有向量与法向量的点积绝对值都小于1e-2，认为点在同一侧
			bool all_small = true;
			for (int k = 0; k < signs.size(); k++)
			{
				if (k == i || k == j)
					continue;
				if (abs(signs(k)) > 1e-2)
				{
					all_small = false;
					break;
				}
			}
			if (all_small)
				return true;

			// 检查所有点积的符号是否相同
			int num1 = 0, num2 = 0;	//正数和负数的个数
			for (int k = 0; k < signs.size(); k++)
			{
				if (k == i || k == j)
					continue;
				if (signs(k) > 0)
					num1++;
				else if (signs(k) < 0)
					num2++;
				//为0就不管
			}
			if (num1 == 0 || num2 == 0)
				return true;
		}
	}
	return false;
}


double pointToSegDist(Eigen::Vector3d& a, Eigen::Vector3d& b, Eigen::Vector3d& p)
{
	// d = |ap| if r <= 0
	//     |bp| if r >= 1
	//     |cp| otherwise
	// r = (ab dot ap)/pow(|ab|,2)

	Eigen::Vector3d ap = p - a;
	Eigen::Vector3d bp = p - b;
	Eigen::Vector3d ab = b - a;

	double r = ab.dot(ap) / pow(ab.norm(), 2);
	if (r <= 0)
		return ap.norm();
	else if (r >= 1)
		return bp.norm();
	else
		return ap.norm() * sqrt(abs(1 - pow((ab.dot(ap) / (ab.norm() * ap.norm())), 2)));

}





//读取到整体模型骨架后,识别1/48子单元中的骨架和关节
//初始1/48四面体的范围:
//V1: p1=(0,0,0) p2=(1,0,0) p3=(1,0,1) p4=(1,1,1) 特征:x>=y, x>=z, y<=z
//V2: p1=(0,0,0) p2=(0,0,1) p3=(1,0,1) p4=(1,1,1) 特征:x>=y, x<=z, y<=z
//V3: p1=(0,0,0) p2=(0,0,1) p3=(0,1,1) p4=(1,1,1) 特征:x<=y, x<=z, y<=z		 
//其他四面体首先对称到xyz都为正数,然后如果y>z就按y=z平面对称上来,此时y<=z
//这时候在三个四面体内,V1和V2关于x=z对称,V2和V3关于x=y对称
//分别识别是否在V3,是否在V2,从而最后对称到V1中

bool inBaseTetrahedral(Eigen::Vector3d v)
{
	if (v(0) < 0 || v(1) < 0 || v(2) < 0)	//初始四面体x,y,z都>0
		return false;
	if (v(1) > v(2))		//初始四面体y<=z
		return false;
	if (v(0) < v(1))		//排除四面体V3
		return false;
	if (v(0) < v(2))		//排除四面体V2
		return false;
	return true;
}

void getPanettaSymmetryPoint(Eigen::Vector3d& v)
{
	double x = abs(v[0]);
	double y = abs(v[1]);
	double z = abs(v[2]);
	if (y > z)
	{
		double temp = y;
		y = z;
		z = temp;
	}
	if (x < y)
	{
		double temp = x;
		x = y;
		y = temp;
	}
	if (x < z)
	{
		double temp = x;
		x = z;
		z = temp;
	}
	v[0] = x, v[1] = y, v[2] = z;
}


bool isEquation(Eigen::Vector3d line1, Eigen::Vector3d line2)
{
	double eps = 1e-5;				//极小值
	if (fabs(line1(0) - line2(0)) < eps && fabs(line1(1) - line2(1)) < eps && fabs(line1(2) - line2(2)) < eps)return true;
	return false;
}

bool Truss::edgesCanMerge(Edge e1, Edge e2, Edge& e, Eigen::Vector3d& e_joint)
{
	Eigen::Vector3d line1 = e1.v.row(1) - e1.v.row(0);
	Eigen::Vector3d line2 = e2.v.row(1) - e2.v.row(0);
	line1.normalize(); 
	line2.normalize(); 
	if (isEquation(line1,line2)||isEquation(line1,-line2))
	{
		if (isEquation(e1.v.row(0), e2.v.row(0)) && e1.r_inner == e2.r_inner && e1.r_outer == e2.r_outer && e1.r_index == e2.r_index)
		{
			e_joint = e1.v.row(0);
			e.r_inner = e1.r_inner;
			e.r_outer = e1.r_outer;
			e.r_index = e1.r_index;
			e.v.resize(2, 3);
			e.v.row(0) = e1.v.row(1);
			e.v.row(1) = e2.v.row(1);
			return true;
		}
		else if (isEquation(e1.v.row(0), e2.v.row(1)) && e1.r_inner == e2.r_inner && e1.r_outer == e2.r_outer && e1.r_index == e2.r_index)
		{
			e_joint = e1.v.row(0);
			e.r_inner = e1.r_inner;
			e.r_outer = e1.r_outer;
			e.r_index = e1.r_index;
			e.v.resize(2, 3);
			e.v.row(0) = e1.v.row(1);
			e.v.row(1) = e2.v.row(0);
			return true;
		}
		else if (isEquation(e1.v.row(1), e2.v.row(0)) && e1.r_inner == e2.r_inner && e1.r_outer == e2.r_outer && e1.r_index == e2.r_index)
		{
			e_joint = e1.v.row(1);
			e.r_inner = e1.r_inner;
			e.r_outer = e1.r_outer;
			e.r_index = e1.r_index;
			e.v.resize(2, 3);
			e.v.row(0) = e1.v.row(0);
			e.v.row(1) = e2.v.row(1);
			return true;
		}
		else if (isEquation(e1.v.row(1), e2.v.row(1)) && e1.r_inner == e2.r_inner && e1.r_outer == e2.r_outer && e1.r_index == e2.r_index)
		{
			e_joint = e1.v.row(1);
			e.r_inner = e1.r_inner;
			e.r_outer = e1.r_outer;
			e.r_index = e1.r_index;
			e.v.resize(2, 3);
			e.v.row(0) = e1.v.row(0);
			e.v.row(1) = e2.v.row(0);
			return true;
		}
	}
	return false;
}


void Truss::getTopology(std::vector<Edge> edge_list)
{
	nodes.clear();
	edges.clear();
	//对edge_list按对称规则去重
	base_edges.clear();
	base_joints.clear();
	for (int i = 0; i < edge_list.size(); i++)
	{
		Eigen::Vector3d v1 = edge_list[i].v.row(0), v2 = edge_list[i].v.row(1);
		//判断是否在初始1/48四面体中
		if (inBaseTetrahedral(v1) && inBaseTetrahedral(v2))
		{
			//加入到base_edges中
			Edge base_e = edge_list[i];
			base_e.r_index = base_edges.size();
			base_edges.push_back(base_e);

			//加入到base_joints中
			bool haveV1 = false;
			for (int j = 0; j < base_joints.size(); j++)
			{
				if (isEquation(v1, base_joints[j].center))
					haveV1 = true;
			}
			if (!haveV1)
			{
				Joint j;
				j.center = v1;
				j.smooth_factor = smooth_factor;
				j.associated_edges.clear();
				j.smf_index = base_joints.size();
				base_joints.push_back(j);
			}
			bool haveV2 = false;
			for (int j = 0; j < base_joints.size(); j++)
			{
				if (isEquation(v2, base_joints[j].center))
					haveV2 = true;
			}
			if (!haveV2)
			{
				Joint j;
				j.center = v2;
				j.smooth_factor = smooth_factor;
				j.associated_edges.clear();
				j.smf_index = base_joints.size();
				base_joints.push_back(j);
			}

		}
	}
	//将所有边的半径映射到base_edges中
	for (int i = 0; i < edge_list.size(); i++)
	{
		Eigen::Vector3d v1 = edge_list[i].v.row(0), v2 = edge_list[i].v.row(1);
		getPanettaSymmetryPoint(v1);
		getPanettaSymmetryPoint(v2);
		for (int j = 0; j < base_edges.size(); j++)
		{
			if (isEquation(v1, base_edges[j].v.row(0)) && isEquation(v2, base_edges[j].v.row(1)))
			{
				edge_list[i].r_index = j;
				break;
			}
			else if (isEquation(v2, base_edges[j].v.row(0)) && isEquation(v1, base_edges[j].v.row(1)))
			{
				edge_list[i].r_index = j;
				break;
			}
		}
		//应该能找到,找不到就报错
		if (edge_list[i].r_index == -1)
		{
			throw std::runtime_error("Skeleton is incompatible with the requested symmetry: Error in getTopology() edges->base_edges!");
		}

	}
	edges = edge_list;

	//将所有点加入到nodes中
	for (int i = 0; i < edge_list.size(); i++)
	{
		Eigen::Vector3d v1 = edge_list[i].v.row(0), v2 = edge_list[i].v.row(1);
		//加入到nodes中
		bool haveV1 = false;
		for (int j = 0; j < nodes.size(); j++)
		{
			if (isEquation(v1, nodes[j].center))
				haveV1 = true;
		}
		if (!haveV1)
		{
			Joint j;
			j.center = v1;
			j.smooth_factor = smooth_factor;
			j.associated_edges.clear();
			j.smf_index = -1;
			nodes.push_back(j);
		}
		bool haveV2 = false;
		for (int j = 0; j < nodes.size(); j++)
		{
			if (isEquation(v2, nodes[j].center))
				haveV2 = true;
		}
		if (!haveV2)
		{
			Joint j;
			j.center = v2;
			j.smooth_factor = smooth_factor;
			j.associated_edges.clear();
			j.smf_index = -1;
			nodes.push_back(j);
		}
	}

	//将所有nodes映射到base_joints
	for (int i = 0; i < nodes.size(); i++)
	{
		Eigen::Vector3d v = nodes[i].center;
		getPanettaSymmetryPoint(v);
		for(int j = 0;j < base_joints.size(); j++)
		{
			if (isEquation(v, base_joints[j].center))
			{
				nodes[i].smf_index = j;
				break;
			}
		}
		//应该能找到, 找不到就报错
		if (nodes[i].smf_index == -1)
		{
			throw std::runtime_error("Skeleton is incompatible with the requested symmetry: Error in getTopology() nodes->base_joints!");
		}
	}

	//计算所有nodes的相邻边(base_joints不用计算相邻边,和nodes中部分点重复了)
	for (int i = 0; i < edge_list.size(); i++)
	{
		Eigen::Vector3d v1 = edge_list[i].v.row(0), v2 = edge_list[i].v.row(1);
		for (int j = 0; j < nodes.size(); j++)
		{
			if (isEquation(v1, nodes[j].center) || isEquation(v2, nodes[j].center))
			{
				nodes[j].associated_edges.push_back(i);
			}
		}
	}
	

	//edges,两条边若能合成一条边,则组合:两条边需要在同一个base_edge
	//这里没考虑三条边ABC连在一起的情况,应该没有这种情况
	//需要考虑这两条边的交点有其他边的情况
	std::vector<Edge> new_edges;

	for (int i = 0; i < edges.size(); i++)
	{
		bool canMerge = false;
		if (edges[i].r_outer == -1)
			continue;
		for (int j = i + 1; j < edges.size(); j++)
		{
			Edge e;
			Eigen::Vector3d e_joint;	//两条边的交点
			if (edgesCanMerge(edges[i], edges[j], e, e_joint))
			{
				//e_joint在nodes中的相邻边为2
				for (int k = 0; k < nodes.size(); k++)
				{
					if (isEquation(nodes[k].center, e_joint) && nodes[k].associated_edges.size()==2)
					{
						canMerge = true;
						break;
					}
				}
				if (canMerge)
				{
					new_edges.push_back(e);
					edges[j].r_outer = -1;
					break;
				}
			}
		}
		if (!canMerge)
			new_edges.push_back(edges[i]);
	}
	edges = new_edges;


	//再次计算所有nodes的相邻边,并去除空nodes
	std::vector<Joint> new_nodes;

	for (int i = 0; i < nodes.size(); i++)
		nodes[i].associated_edges.clear();

	for (int i = 0; i < edges.size(); i++)
	{
		Eigen::Vector3d v1 = edges[i].v.row(0), v2 = edges[i].v.row(1);
		for (int j = 0; j < nodes.size(); j++)
		{
			if (isEquation(v1, nodes[j].center) || isEquation(v2, nodes[j].center))
			{
				nodes[j].associated_edges.push_back(i);
			}
		}
	}

	for (int i = 0; i < nodes.size(); i++)
	{
		if (nodes[i].associated_edges.size() > 0)
			new_nodes.push_back(nodes[i]);
	}
	nodes = new_nodes;

	//根据nodes映射到base_joints中
	std::vector<Joint> new_base_joints;
	for (int i = 0; i < base_joints.size(); i++)
	{
		bool needDelete = true;
		for (int j = 0; j < nodes.size(); j++)
		{
			if (nodes[j].smf_index == i)
			{
				needDelete = false;
				break;
			}
		}
		if (!needDelete)
		{
			//如果base_joint是原位置,直接插入
			if (new_base_joints.size() == i)
			{
				new_base_joints.push_back(base_joints[i]);
			}
			//否则还要修改nodes中的smf_index索引
			else
			{
				int new_smf_index = new_base_joints.size();
				new_base_joints.push_back(base_joints[i]);
				for (int j = 0; j < nodes.size(); j++)
				{
					if (nodes[j].smf_index == i)
						nodes[j].smf_index = new_smf_index;
				}
			}
			
		}
	}
	base_joints = new_base_joints;
	//对所有nodes计算相邻边是否在空间同一侧
	//如果对应base_joints有赋值,则直接取即可
	std::vector<bool> same_side_already;
	for (int i = 0; i < base_joints.size(); i++)
		same_side_already.push_back(false);
	for (int i = 0; i < nodes.size(); i++)
	{
		//首先查找对应base_joints
		int base_index = nodes[i].smf_index;
		if (same_side_already[base_index])
		{
			nodes[i].is_same_side = base_joints[base_index].is_same_side;
			continue;
		}
		//否则计算并赋值给对应base_joints
		Eigen::Vector3d joint_center = nodes[i].center;
		int points_size = nodes[i].associated_edges.size(); 
		Eigen::MatrixXd points(points_size, 3);
		for (int j = 0; j < points_size; j++)
		{
			Edge e = edges[nodes[i].associated_edges[j]];
			if (isEquation(joint_center, e.v.row(0)))
				points.row(j) = e.v.row(1);
			else if (isEquation(joint_center, e.v.row(1)))
				points.row(j) = e.v.row(0);
			else
				std::cout<< "Error in getTopology() calculate is_same_side" << std::endl;
		}
		nodes[i].is_same_side = edgesAreSameSide(joint_center, points);
		base_joints[base_index].is_same_side = nodes[i].is_same_side;
		same_side_already[base_index] = true;
	}
}


void Truss::getSkeletonFromPanetta(std::string path)
{
	Eigen::MatrixXd csv;
	if (!elib::readCSV(path, csv) || csv.cols() != 3 || csv.rows() < 2 || csv.rows() % 2 != 0)
        throw std::runtime_error("Skeleton must contain numeric CSV rows in endpoint pairs (3 columns): " + path);
    if (!csv.allFinite() || csv.cwiseAbs().maxCoeff() > 1.0 + 1e-8)
        throw std::runtime_error("Skeleton coordinates must be finite and normalized to [-1, 1]: " + path);
    for (Eigen::Index i = 0; i < csv.rows(); i += 2)
        if ((csv.row(i) - csv.row(i + 1)).norm() < 1e-10)
            throw std::runtime_error("Skeleton contains a zero-length edge: " + path);

	std::vector<Edge> edgelist;
	for (int i = 0; i < csv.rows(); i += 2)
	{
		Eigen::Vector3d p1(csv(i, 0), csv(i, 1), csv(i, 2));
		Eigen::Vector3d p2(csv(i + 1, 0), csv(i + 1, 1), csv(i + 1, 2));

		Edge e1;
		e1.v.resize(2, 3);
		e1.v.row(0) = p1 * (cell_length / 2);
		e1.v.row(1) = p2 * (cell_length / 2);
		e1.r_outer = r_outer;
		e1.r_inner = r_inner;
		e1.r_index = -1;
		edgelist.push_back(e1);
	}

	getTopology(edgelist);

	//这里省掉了PSL中的控制点和厚度thickness的计算
}


bool inBaseCube(Eigen::Vector3d v)
{
	if (v(0) >= 0 && v(1) >= 0 & v(2) >= 0)
		return true;
	return false;
}

void getETHSymmetryPoint(Eigen::Vector3d& v)
{
	v[0] = abs(v[0]);
	v[1] = abs(v[1]);
	v[2] = abs(v[2]);
}

void Truss::getTopologyETH(std::vector<Edge> edge_list)
{
	nodes.clear();
	edges.clear();
	//对edge_list按对称规则去重
	base_edges.clear();
	base_joints.clear();
	for (int i = 0; i < edge_list.size(); i++)
	{
		Eigen::Vector3d v1 = edge_list[i].v.row(0), v2 = edge_list[i].v.row(1);
		//判断是否在初始1/8立方体中
		if (inBaseCube(v1) && inBaseCube(v2))
		{
			//加入到base_edges中
			Edge base_e = edge_list[i];
			base_e.r_index = base_edges.size();
			base_edges.push_back(base_e);

			//加入到base_joints中
			bool haveV1 = false;
			for (int j = 0; j < base_joints.size(); j++)
			{
				if (isEquation(v1, base_joints[j].center))
					haveV1 = true;
			}
			if (!haveV1)
			{
				Joint j;
				j.center = v1;
				j.smooth_factor = smooth_factor;
				j.associated_edges.clear();
				j.smf_index = base_joints.size();
				base_joints.push_back(j);
			}
			bool haveV2 = false;
			for (int j = 0; j < base_joints.size(); j++)
			{
				if (isEquation(v2, base_joints[j].center))
					haveV2 = true;
			}
			if (!haveV2)
			{
				Joint j;
				j.center = v2;
				j.smooth_factor = smooth_factor;
				j.associated_edges.clear();
				j.smf_index = base_joints.size();
				base_joints.push_back(j);
			}

		}
	}
	//将所有边的半径映射到base_edges中
	for (int i = 0; i < edge_list.size(); i++)
	{
		Eigen::Vector3d v1 = edge_list[i].v.row(0), v2 = edge_list[i].v.row(1);
		getETHSymmetryPoint(v1);
		getETHSymmetryPoint(v2);
		for (int j = 0; j < base_edges.size(); j++)
		{
			if (isEquation(v1, base_edges[j].v.row(0)) && isEquation(v2, base_edges[j].v.row(1)))
			{
				edge_list[i].r_index = j;
				break;
			}
			else if (isEquation(v2, base_edges[j].v.row(0)) && isEquation(v1, base_edges[j].v.row(1)))
			{
				edge_list[i].r_index = j;
				break;
			}
		}
		//应该能找到,找不到就报错
		if (edge_list[i].r_index == -1)
		{
			throw std::runtime_error("Skeleton is incompatible with the requested symmetry: Error in getTopology() edges->base_edges!");
		}

	}
	edges = edge_list;

	//将所有点加入到nodes中
	for (int i = 0; i < edge_list.size(); i++)
	{
		Eigen::Vector3d v1 = edge_list[i].v.row(0), v2 = edge_list[i].v.row(1);
		//加入到nodes中
		bool haveV1 = false;
		for (int j = 0; j < nodes.size(); j++)
		{
			if (isEquation(v1, nodes[j].center))
				haveV1 = true;
		}
		if (!haveV1)
		{
			Joint j;
			j.center = v1;
			j.smooth_factor = smooth_factor;
			j.associated_edges.clear();
			j.smf_index = -1;
			nodes.push_back(j);
		}
		bool haveV2 = false;
		for (int j = 0; j < nodes.size(); j++)
		{
			if (isEquation(v2, nodes[j].center))
				haveV2 = true;
		}
		if (!haveV2)
		{
			Joint j;
			j.center = v2;
			j.smooth_factor = smooth_factor;
			j.associated_edges.clear();
			j.smf_index = -1;
			nodes.push_back(j);
		}
	}

	//将所有nodes映射到base_joints
	for (int i = 0; i < nodes.size(); i++)
	{
		Eigen::Vector3d v = nodes[i].center;
		getETHSymmetryPoint(v);
		for (int j = 0; j < base_joints.size(); j++)
		{
			if (isEquation(v, base_joints[j].center))
			{
				nodes[i].smf_index = j;
				break;
			}
		}
		//应该能找到, 找不到就报错
		if (nodes[i].smf_index == -1)
		{
			throw std::runtime_error("Skeleton is incompatible with the requested symmetry: Error in getTopology() nodes->base_joints!");
		}
	}

	//计算所有nodes的相邻边(base_joints不用计算相邻边,和nodes中部分点重复了)
	for (int i = 0; i < edge_list.size(); i++)
	{
		Eigen::Vector3d v1 = edge_list[i].v.row(0), v2 = edge_list[i].v.row(1);
		for (int j = 0; j < nodes.size(); j++)
		{
			if (isEquation(v1, nodes[j].center) || isEquation(v2, nodes[j].center))
			{
				nodes[j].associated_edges.push_back(i);
			}
		}
	}


	//edges,两条边若能合成一条边,则组合:两条边需要在同一个base_edge
	//这里没考虑三条边ABC连在一起的情况,应该没有这种情况
	//需要考虑这两条边的交点有其他边的情况
	std::vector<Edge> new_edges;

	for (int i = 0; i < edges.size(); i++)
	{
		bool canMerge = false;
		if (edges[i].r_outer == -1)
			continue;
		for (int j = i + 1; j < edges.size(); j++)
		{
			Edge e;
			Eigen::Vector3d e_joint;	//两条边的交点
			if (edgesCanMerge(edges[i], edges[j], e, e_joint))
			{
				//e_joint在nodes中的相邻边为2
				for (int k = 0; k < nodes.size(); k++)
				{
					if (isEquation(nodes[k].center, e_joint) && nodes[k].associated_edges.size() == 2)
					{
						canMerge = true;
						break;
					}
				}
				if (canMerge)
				{
					new_edges.push_back(e);
					edges[j].r_outer = -1;
					break;
				}
			}
		}
		if (!canMerge)
			new_edges.push_back(edges[i]);
	}
	edges = new_edges;


	//再次计算所有nodes的相邻边,并去除空nodes
	std::vector<Joint> new_nodes;

	for (int i = 0; i < nodes.size(); i++)
		nodes[i].associated_edges.clear();

	for (int i = 0; i < edges.size(); i++)
	{
		Eigen::Vector3d v1 = edges[i].v.row(0), v2 = edges[i].v.row(1);
		for (int j = 0; j < nodes.size(); j++)
		{
			if (isEquation(v1, nodes[j].center) || isEquation(v2, nodes[j].center))
			{
				nodes[j].associated_edges.push_back(i);
			}
		}
	}

	for (int i = 0; i < nodes.size(); i++)
	{
		if (nodes[i].associated_edges.size() > 0)
			new_nodes.push_back(nodes[i]);
	}
	nodes = new_nodes;

	//根据nodes映射到base_joints中
	std::vector<Joint> new_base_joints;
	for (int i = 0; i < base_joints.size(); i++)
	{
		bool needDelete = true;
		for (int j = 0; j < nodes.size(); j++)
		{
			if (nodes[j].smf_index == i)
			{
				needDelete = false;
				break;
			}
		}
		if (!needDelete)
		{
			//如果base_joint是原位置,直接插入
			if (new_base_joints.size() == i)
			{
				new_base_joints.push_back(base_joints[i]);
			}
			//否则还要修改nodes中的smf_index索引
			else
			{
				int new_smf_index = new_base_joints.size();
				new_base_joints.push_back(base_joints[i]);
				for (int j = 0; j < nodes.size(); j++)
				{
					if (nodes[j].smf_index == i)
						nodes[j].smf_index = new_smf_index;
				}
			}

		}
	}
	base_joints = new_base_joints;
	//对所有nodes计算相邻边是否在空间同一侧
	for (int i = 0; i < nodes.size(); i++)
	{
		Eigen::Vector3d joint_center = nodes[i].center;
		int points_size = nodes[i].associated_edges.size();
		Eigen::MatrixXd points(points_size, 3);
		for (int j = 0; j < points_size; j++)
		{
			Edge e = edges[nodes[i].associated_edges[j]];
			if (isEquation(joint_center, e.v.row(0)))
				points.row(j) = e.v.row(1);
			else if (isEquation(joint_center, e.v.row(1)))
				points.row(j) = e.v.row(0);
			else
				std::cout << "Error in getTopology() calculate is_same_side" << std::endl;
		}
		nodes[i].is_same_side = edgesAreSameSide(joint_center, points);
	}
}

void Truss::getSkeletonFromETH(std::string path)
{
	Eigen::MatrixXd csv;
	if (!elib::readCSV(path, csv) || csv.cols() != 3 || csv.rows() < 2 || csv.rows() % 2 != 0)
        throw std::runtime_error("Skeleton must contain numeric CSV rows in endpoint pairs (3 columns): " + path);
    if (!csv.allFinite() || csv.cwiseAbs().maxCoeff() > 1.0 + 1e-8)
        throw std::runtime_error("Skeleton coordinates must be finite and normalized to [-1, 1]: " + path);
    for (Eigen::Index i = 0; i < csv.rows(); i += 2)
        if ((csv.row(i) - csv.row(i + 1)).norm() < 1e-10)
            throw std::runtime_error("Skeleton contains a zero-length edge: " + path);

	std::vector<Edge> edgelist;
	for (int i = 0; i < csv.rows(); i += 2)
	{
		Eigen::Vector3d p1(csv(i, 0), csv(i, 1), csv(i, 2));
		Eigen::Vector3d p2(csv(i + 1, 0), csv(i + 1, 1), csv(i + 1, 2));

		Edge e1;
		e1.v.resize(2, 3);
		e1.v.row(0) = p1 * (cell_length / 2);
		e1.v.row(1) = p2 * (cell_length / 2);
		e1.r_outer = r_outer;
		e1.r_inner = r_inner;
		e1.r_index = -1;
		edgelist.push_back(e1);
	}

	getTopologyETH(edgelist);

	//这里省掉了PSL中的控制点和厚度thickness的计算
}

bool Truss::setSolidVolumeFraction(double vof)
{
	//计算当前体积分数
	genSolidTrussVoxel();
	double vof_1 = (double)voxel.sum() / (double)voxel.rows();
	//std::cout << "初始体积分数为" << vof_1 << ",半径为" << r_outer << std::endl;

	//将d_outer加一点点,然后计算斜率
	setROuter(r_outer + 0.5, -1);
	genSolidTrussVoxel();
	double vof_2 = (double)voxel.sum() / (double)voxel.rows();
	//std::cout << "初始体积分数为" << vof_2 << ",半径为" << r_outer << std::endl;

	double k = ((vof_2 - vof_1) * 1.0) / 0.5;//斜率
	if (!std::isfinite(k)) return false;
	double ini_r = std::abs(k) < 1e-14 ? r_outer : ((vof - vof_2) / k + r_outer);//粗略计算目标相对密度的半径
	if (ini_r < 0)
		ini_r = 1.0;		//从最小开始
	
	//迭代求解精确的半径,确保体积分数在+-0.5%
	setROuter(ini_r, -1);
	genSolidTrussVoxel();
	double vof_current = (double)voxel.sum() / (double)voxel.rows();
	//std::cout << "迭代开始时体积分数为" << vof_current << ",半径为" << ini_r << std::endl;

	double factor = 1.0;	//步进长度因子
	int fit_iterations = 0;
	while (std::abs(vof - vof_current) > 0.005)
	{
		if (factor < 1e-10 || ++fit_iterations > 1000) return false;
		if (vof_current < vof)
		{
			//一点点增加
			while (vof - vof_current > 0.005)
			{
				if (++fit_iterations > 1000) return false;
				ini_r += 0.5 * factor;
				setROuter(ini_r, -1);
				genSolidTrussVoxel();
				vof_current = (double)voxel.sum() / (double)voxel.rows();
			}
		}
		else if (vof_current > vof)
		{
			//一点点减少
			while (vof - vof_current < -0.005)
			{
				if (++fit_iterations > 1000) return false;
				ini_r -= 0.5 * factor;
				setROuter(ini_r, -1);
				genSolidTrussVoxel();
				vof_current = (double)voxel.sum() / (double)voxel.rows();
			}
		}
		else {}
		factor = factor * 0.1;
		//std::cout << "迭代时体积分数为" << vof_current << ",半径为" << ini_r << std::endl;
	}
	std::cout << "最终体积分数为" << vof_current << ",半径为" << ini_r << std::endl;
	if (ini_r < 0.8)
		return false;
	return true;
}

void Truss::genSolidTrussVoxel()
{
	double cell_length_half = cell_length / 2;	//模型大小的一半,与一些数值参数有关

	//重建范围
	Eigen::RowVector3d VminCut = { -1 * cell_length_half, -1 * cell_length_half, -1 * cell_length_half };
	Eigen::RowVector3d VmaxCut = { cell_length_half, cell_length_half, cell_length_half};

	Eigen::RowVector3i res = (resolution * ((VmaxCut - VminCut) / (VmaxCut - VminCut).maxCoeff())).cast<int>();

	grid.resize(res(0) * res(1) * res(2), 3);

	//构建采样点的坐标
	const auto lerpCut = [&](const int di, const int d)->double
	{return VminCut(d) + (double)di / (double)(res(d) - 1) * (VmaxCut(d) - VminCut(d)); };

#pragma omp parallel for
	for (int zi = 0; zi < res(2); zi++)
	{
		const double zc = lerpCut(zi, 2);
		for (int yi = 0; yi < res(1); yi++)
		{
			const double yc = lerpCut(yi, 1);
			for (int xi = 0; xi < res(0); xi++)
			{
				const double xc = lerpCut(xi, 0);
				grid.row(xi + res(0) * (yi + res(1) * zi)) = Eigen::RowVector3d(xc, yc, zc);
			}
		}
	}

	//std::vector<Edge> edge_list = base_edges;
	std::vector<Edge> edge_list = edges;
	std::vector<Joint> joint_list = nodes;
	//double smf = smooth_factor;

	//计算距离场的公式
	//看PPL的计算
	//将采样点xyz映射到子单元1/8或者1/48中,然后计算距离,只需要关照子单元中的几条杆
	const auto& function = [edge_list, joint_list, cell_length_half](double x, double y, double z)->double {
		//超出边界,赋一个比较大的正值
		if (abs(x) > cell_length_half || abs(y) > cell_length_half || abs(z) > cell_length_half)
		{
			return 100;
		}
		//将采样点映射到子单元
		Eigen::Vector3d p = { x,y,z };	
		getPanettaSymmetryPoint(p);

		std::vector<double> result;	//存放每个joint的KS函数值,到最后再套一层总的KS函数

		for (int i = 0; i < joint_list.size(); i++)
		{
			double result_joint = 0;		//当前joint的KS函数值
			double smf = joint_list[i].smooth_factor;
			for (int j = 0; j < joint_list[i].associated_edges.size(); j++)
			{
				int index = joint_list[i].associated_edges[j];	//相邻边在edge_list中的索引
				Eigen::Vector3d v1 = edge_list[index].v.row(0), v2 = edge_list[index].v.row(1);
				double dis = pointToSegDist(v1, v2, p) - edge_list[index].r_outer;
				result_joint += exp(-1 / smf * dis);	//smf是公式中的1/s_j
			}
			result_joint = -smf * log(result_joint);
			result.push_back(result_joint);
		}

		double sum_KS = 0;	//总的KS函数值
		double sum_smf = 0.5;
		for (int i = 0; i < result.size(); i++)
		{
			sum_KS += exp(-1 / sum_smf * result[i]);
		}
		sum_KS = -sum_smf * log(sum_KS);
		return sum_KS;
	};
	voxel.resize(grid.rows());
	voxel.setZero();
#pragma omp parallel for
	for (int i = 0; i < grid.rows(); i++) {
		double indicator = function(grid(i, 0), grid(i, 1), grid(i, 2));
		if (indicator <= 0)
		{
			voxel(i) = 1;
		}
		else
			voxel(i) = 0;
	}
}

void Truss::genSolidTrussMesh()
{
	double cell_length_half = cell_length / 2;	//模型大小的一半,与一些数值参数有关

	//重建范围
	Eigen::RowVector3d Vmin = { -1 * cell_length_half - 0.1, -1 * cell_length_half - 0.1, -1 * cell_length_half - 0.1 };
	Eigen::RowVector3d Vmax = { cell_length_half + 0.1, cell_length_half + 0.1, cell_length_half + 0.1 };

	Eigen::RowVector3i res_mesh = (resolution_mesh * ((Vmax - Vmin) / (Vmax - Vmin).maxCoeff())).cast<int>();

	Eigen::MatrixXd GV(res_mesh(0) * res_mesh(1) * res_mesh(2), 3);

	//构建采样点的坐标
	const auto lerp = [&](const int di, const int d)->double
	{return Vmin(d) + (double)di / (double)(res_mesh(d) - 1) * (Vmax(d) - Vmin(d)); };

#pragma omp parallel for
	for (int zi = 0; zi < res_mesh(2); zi++)
	{
		const double z = lerp(zi, 2);
		for (int yi = 0; yi < res_mesh(1); yi++)
		{
			const double y = lerp(yi, 1);
			for (int xi = 0; xi < res_mesh(0); xi++)
			{
				const double x = lerp(xi, 0);
				GV.row(xi + res_mesh(0) * (yi + res_mesh(1) * zi)) = Eigen::RowVector3d(x, y, z);
			}
		}
	}

	//std::vector<Edge> edge_list = base_edges;
	std::vector<Edge> edge_list = edges;
	std::vector<Joint> joint_list = nodes;
	//double smf = smooth_factor;

	//计算距离场的公式
	//看PPL的计算
	//将采样点xyz映射到子单元1/8或者1/48中,然后计算距离,只需要关照子单元中的几条杆
	const auto& function = [edge_list, joint_list, cell_length_half](double x, double y, double z)->double {
		//超出边界,赋一个比较大的正值
		if (abs(x) > cell_length_half || abs(y) > cell_length_half || abs(z) > cell_length_half)
		{
			return 100;
		}
		//将采样点映射到子单元
		Eigen::Vector3d p = { x,y,z };
		getPanettaSymmetryPoint(p);

		std::vector<double> result;	//存放每个joint的KS函数值,到最后再套一层总的KS函数

		for (int i = 0; i < joint_list.size(); i++)
		{
			double result_joint = 0;		//当前joint的KS函数值
			double smf = joint_list[i].smooth_factor;
			for (int j = 0; j < joint_list[i].associated_edges.size(); j++)
			{
				int index = joint_list[i].associated_edges[j];	//相邻边在edge_list中的索引
				Eigen::Vector3d v1 = edge_list[index].v.row(0), v2 = edge_list[index].v.row(1);
				double dis = pointToSegDist(v1, v2, p) - edge_list[index].r_outer;
				result_joint += exp(-1 / smf * dis);	//smf是公式中的1/s_j
			}
			result_joint = -smf * log(result_joint);
			result.push_back(result_joint);
		}

		double sum_KS = 0;	//总的KS函数值
		double sum_smf = 0.5;
		for (int i = 0; i < result.size(); i++)
		{
			sum_KS += exp(-1 / sum_smf * result[i]);
		}
		sum_KS = -sum_smf * log(sum_KS);
		return sum_KS;
	};

	Eigen::VectorXd B(GV.rows());
#pragma omp parallel for
	for (int i = 0; i < GV.rows(); i++) {
		double indicator = function(GV(i, 0), GV(i, 1), GV(i, 2));
		B(i) = indicator;
	}
	igl::copyleft::marching_cubes(B, GV, res_mesh(0), res_mesh(1), res_mesh(2), V, F);

	//平移中心到原点
	if (V.rows() == 0 || !V.allFinite())
		throw std::runtime_error("No finite isosurface; check the full-cell skeleton, radius and mesh resolution");
	Eigen::Vector3d minV = V.colwise().minCoeff();
	Eigen::Vector3d maxV = V.colwise().maxCoeff();

	Eigen::Vector3d move = minV + (maxV - minV) / 2.0;
	for (int j = 0; j < V.rows(); j++)
	{
		V(j, 0) -= move(0);
		V(j, 1) -= move(1);
		V(j, 2) -= move(2);
	}

	//缩放为[-1,-1,-1]到[1,1,1]的边长为2的立方体
	Eigen::Vector3d factor = (maxV - minV) / 2.0;
	double xfac = (1.0 / factor(0));
	double yfac = (1.0 / factor(1));
	double zfac = (1.0 / factor(2));
	for (int j = 0; j < V.rows(); j++)
	{

		V(j, 0) *= xfac;
		V(j, 1) *= yfac;
		V(j, 2) *= zfac;
	}
}

bool Truss::setSolidVolumeFractionETH(double vof)
{
	//计算当前体积分数
	genSolidTrussVoxelETH();
	double vof_1 = (double)voxel.sum() / (double)voxel.rows();
	//std::cout << "初始体积分数为" << vof_1 << ",半径为" << r_outer << std::endl;

	//将d_outer加一点点,然后计算斜率
	setROuter(r_outer + 0.5, -1);
	genSolidTrussVoxelETH();
	double vof_2 = (double)voxel.sum() / (double)voxel.rows();
	//std::cout << "初始体积分数为" << vof_2 << ",半径为" << r_outer << std::endl;

	double k = ((vof_2 - vof_1) * 1.0) / 0.5;//斜率
	if (!std::isfinite(k)) return false;
	double ini_r = std::abs(k) < 1e-14 ? r_outer : ((vof - vof_2) / k + r_outer);//粗略计算目标相对密度的半径
	if (ini_r < 0)
		ini_r = 1.0;		//从最小开始

	//迭代求解精确的半径,确保体积分数在+-0.5%
	setROuter(ini_r, -1);
	genSolidTrussVoxelETH();
	double vof_current = (double)voxel.sum() / (double)voxel.rows();
	//std::cout << "迭代开始时体积分数为" << vof_current << ",半径为" << ini_r << std::endl;

	double factor = 1.0;	//步进长度因子
	int fit_iterations = 0;
	while (std::abs(vof - vof_current) > 0.005)
	{
		if (factor < 1e-10 || ++fit_iterations > 1000) return false;
		if (vof_current < vof)
		{
			//一点点增加
			while (vof - vof_current > 0.005)
			{
				if (++fit_iterations > 1000) return false;
				ini_r += 0.5 * factor;
				setROuter(ini_r, -1);
				genSolidTrussVoxelETH();
				vof_current = (double)voxel.sum() / (double)voxel.rows();
			}
		}
		else if (vof_current > vof)
		{
			//一点点减少
			while (vof - vof_current < -0.005)
			{
				if (++fit_iterations > 1000) return false;
				ini_r -= 0.5 * factor;
				setROuter(ini_r, -1);
				genSolidTrussVoxelETH();
				vof_current = (double)voxel.sum() / (double)voxel.rows();
			}
		}
		else {}
		factor = factor * 0.1;
		//std::cout << "迭代时体积分数为" << vof_current << ",半径为" << ini_r << std::endl;
	}
	std::cout << "最终体积分数为" << vof_current << ",半径为" << ini_r << std::endl;
	if (ini_r < 0.8)
		return false;
	return true;
}

void Truss::genSolidTrussVoxelETH()
{
	double cell_length_half = cell_length / 2;	//模型大小的一半,与一些数值参数有关

	//重建范围
	Eigen::RowVector3d VminCut = { -1 * cell_length_half, -1 * cell_length_half, -1 * cell_length_half };
	Eigen::RowVector3d VmaxCut = { cell_length_half, cell_length_half, cell_length_half };

	Eigen::RowVector3i res = (resolution * ((VmaxCut - VminCut) / (VmaxCut - VminCut).maxCoeff())).cast<int>();

	grid.resize(res(0) * res(1) * res(2), 3);

	//构建采样点的坐标
	const auto lerpCut = [&](const int di, const int d)->double
	{return VminCut(d) + (double)di / (double)(res(d) - 1) * (VmaxCut(d) - VminCut(d)); };

#pragma omp parallel for
	for (int zi = 0; zi < res(2); zi++)
	{
		const double zc = lerpCut(zi, 2);
		for (int yi = 0; yi < res(1); yi++)
		{
			const double yc = lerpCut(yi, 1);
			for (int xi = 0; xi < res(0); xi++)
			{
				const double xc = lerpCut(xi, 0);
				grid.row(xi + res(0) * (yi + res(1) * zi)) = Eigen::RowVector3d(xc, yc, zc);
			}
		}
	}

	//std::vector<Edge> edge_list = base_edges;
	std::vector<Edge> edge_list = edges;
	std::vector<Joint> joint_list = nodes;
	//double smf = smooth_factor;

	//计算距离场的公式
	//看PPL的计算
	//将采样点xyz映射到子单元1/8或者1/48中,然后计算距离,只需要关照子单元中的几条杆
	const auto& function = [edge_list, joint_list, cell_length_half](double x, double y, double z)->double {
		//超出边界,赋一个比较大的正值
		if (abs(x) > cell_length_half || abs(y) > cell_length_half || abs(z) > cell_length_half)
		{
			return 100;
		}
		//将采样点映射到子单元
		Eigen::Vector3d p = { x,y,z };
		getETHSymmetryPoint(p);

		std::vector<double> result;	//存放每个joint的KS函数值,到最后再套一层总的KS函数

		for (int i = 0; i < joint_list.size(); i++)
		{
			double result_joint = 0;		//当前joint的KS函数值
			double smf = joint_list[i].smooth_factor;
			for (int j = 0; j < joint_list[i].associated_edges.size(); j++)
			{
				int index = joint_list[i].associated_edges[j];	//相邻边在edge_list中的索引
				Eigen::Vector3d v1 = edge_list[index].v.row(0), v2 = edge_list[index].v.row(1);
				double dis = pointToSegDist(v1, v2, p) - edge_list[index].r_outer;
				result_joint += exp(-1 / smf * dis);	//smf是公式中的1/s_j
			}
			result_joint = -smf * log(result_joint);
			result.push_back(result_joint);
		}

		double sum_KS = 0;	//总的KS函数值
		double sum_smf = 0.5;
		for (int i = 0; i < result.size(); i++)
		{
			sum_KS += exp(-1 / sum_smf * result[i]);
		}
		sum_KS = -sum_smf * log(sum_KS);
		return sum_KS;
	};
	voxel.resize(grid.rows());
	voxel.setZero();
#pragma omp parallel for
	for (int i = 0; i < grid.rows(); i++) {
		double indicator = function(grid(i, 0), grid(i, 1), grid(i, 2));
		if (indicator <= 0)
		{
			voxel(i) = 1;
		}
		else
			voxel(i) = 0;
	}
}

void Truss::genSolidTrussMeshETH()
{
	double cell_length_half = cell_length / 2;	//模型大小的一半,与一些数值参数有关

	//重建范围
	Eigen::RowVector3d Vmin = { -1 * cell_length_half - 0.1, -1 * cell_length_half - 0.1, -1 * cell_length_half - 0.1 };
	Eigen::RowVector3d Vmax = { cell_length_half + 0.1, cell_length_half + 0.1, cell_length_half + 0.1 };

	Eigen::RowVector3i res_mesh = (resolution_mesh * ((Vmax - Vmin) / (Vmax - Vmin).maxCoeff())).cast<int>();

	Eigen::MatrixXd GV(res_mesh(0) * res_mesh(1) * res_mesh(2), 3);

	//构建采样点的坐标
	const auto lerp = [&](const int di, const int d)->double
	{return Vmin(d) + (double)di / (double)(res_mesh(d) - 1) * (Vmax(d) - Vmin(d)); };

#pragma omp parallel for
	for (int zi = 0; zi < res_mesh(2); zi++)
	{
		const double z = lerp(zi, 2);
		for (int yi = 0; yi < res_mesh(1); yi++)
		{
			const double y = lerp(yi, 1);
			for (int xi = 0; xi < res_mesh(0); xi++)
			{
				const double x = lerp(xi, 0);
				GV.row(xi + res_mesh(0) * (yi + res_mesh(1) * zi)) = Eigen::RowVector3d(x, y, z);
			}
		}
	}

	//std::vector<Edge> edge_list = base_edges;
	std::vector<Edge> edge_list = edges;
	std::vector<Joint> joint_list = nodes;
	//double smf = smooth_factor;

	//计算距离场的公式
	//看PPL的计算
	//将采样点xyz映射到子单元1/8或者1/48中,然后计算距离,只需要关照子单元中的几条杆
	const auto& function = [edge_list, joint_list, cell_length_half](double x, double y, double z)->double {
		//超出边界,赋一个比较大的正值
		if (abs(x) > cell_length_half || abs(y) > cell_length_half || abs(z) > cell_length_half)
		{
			return 100;
		}
		//将采样点映射到子单元
		Eigen::Vector3d p = { x,y,z };
		getETHSymmetryPoint(p);

		std::vector<double> result;	//存放每个joint的KS函数值,到最后再套一层总的KS函数

		for (int i = 0; i < joint_list.size(); i++)
		{
			double result_joint = 0;		//当前joint的KS函数值
			double smf = joint_list[i].smooth_factor;
			for (int j = 0; j < joint_list[i].associated_edges.size(); j++)
			{
				int index = joint_list[i].associated_edges[j];	//相邻边在edge_list中的索引
				Eigen::Vector3d v1 = edge_list[index].v.row(0), v2 = edge_list[index].v.row(1);
				double dis = pointToSegDist(v1, v2, p) - edge_list[index].r_outer;
				result_joint += exp(-1 / smf * dis);	//smf是公式中的1/s_j
			}
			result_joint = -smf * log(result_joint);
			result.push_back(result_joint);
		}

		double sum_KS = 0;	//总的KS函数值
		double sum_smf = 0.5;
		for (int i = 0; i < result.size(); i++)
		{
			sum_KS += exp(-1 / sum_smf * result[i]);
		}
		sum_KS = -sum_smf * log(sum_KS);
		return sum_KS;
	};

	Eigen::VectorXd B(GV.rows());
#pragma omp parallel for
	for (int i = 0; i < GV.rows(); i++) {
		double indicator = function(GV(i, 0), GV(i, 1), GV(i, 2));
		B(i) = indicator;
	}
	igl::copyleft::marching_cubes(B, GV, res_mesh(0), res_mesh(1), res_mesh(2), V, F);

	//平移中心到原点
	if (V.rows() == 0 || !V.allFinite())
		throw std::runtime_error("No finite isosurface; check the full-cell skeleton, radius and mesh resolution");
	Eigen::Vector3d minV = V.colwise().minCoeff();
	Eigen::Vector3d maxV = V.colwise().maxCoeff();

	Eigen::Vector3d move = minV + (maxV - minV) / 2.0;
	for (int j = 0; j < V.rows(); j++)
	{
		V(j, 0) -= move(0);
		V(j, 1) -= move(1);
		V(j, 2) -= move(2);
	}

	//缩放为[-1,-1,-1]到[1,1,1]的边长为2的立方体
	Eigen::Vector3d factor = (maxV - minV) / 2.0;
	double xfac = (1.0 / factor(0));
	double yfac = (1.0 / factor(1));
	double zfac = (1.0 / factor(2));
	for (int j = 0; j < V.rows(); j++)
	{

		V(j, 0) *= xfac;
		V(j, 1) *= yfac;
		V(j, 2) *= zfac;
	}
}


bool Truss::setNormalVolumeFractionETH(double vof)
{
	//计算当前体积分数
	genNormalTrussVoxelETH();
	double vof_1 = (double)voxel.sum() / (double)voxel.rows();
	//std::cout << "初始体积分数为" << vof_1 << ",半径为" << r_outer << std::endl;

	//将d_outer加一点点,然后计算斜率
	setROuter(r_outer + 0.5, -1);
	genNormalTrussVoxelETH();
	double vof_2 = (double)voxel.sum() / (double)voxel.rows();
	//std::cout << "初始体积分数为" << vof_2 << ",半径为" << r_outer << std::endl;

	double k = ((vof_2 - vof_1) * 1.0) / 0.5;//斜率
	if (!std::isfinite(k)) return false;
	double ini_r = std::abs(k) < 1e-14 ? r_outer : ((vof - vof_2) / k + r_outer);//粗略计算目标相对密度的半径
	if (ini_r < 0)
		ini_r = 1.0;		//从最小开始

	//迭代求解精确的半径,确保体积分数在+-0.5%
	setROuter(ini_r, -1);
	genNormalTrussVoxelETH();
	double vof_current = (double)voxel.sum() / (double)voxel.rows();
	//std::cout << "迭代开始时体积分数为" << vof_current << ",半径为" << ini_r << std::endl;

	double factor = 1.0;	//步进长度因子
	int fit_iterations = 0;
	while (std::abs(vof - vof_current) > 0.005)
	{
		if (factor < 1e-10 || ++fit_iterations > 1000) return false;
		if (vof_current < vof)
		{
			//一点点增加
			while (vof - vof_current > 0.005)
			{
				if (++fit_iterations > 1000) return false;
				ini_r += 0.5 * factor;
				setROuter(ini_r, -1);
				genNormalTrussVoxelETH();
				vof_current = (double)voxel.sum() / (double)voxel.rows();
			}
		}
		else if (vof_current > vof)
		{
			//一点点减少
			while (vof - vof_current < -0.005)
			{
				if (++fit_iterations > 1000) return false;
				ini_r -= 0.5 * factor;
				setROuter(ini_r, -1);
				genNormalTrussVoxelETH();
				vof_current = (double)voxel.sum() / (double)voxel.rows();
			}
		}
		else {}
		factor = factor * 0.1;
		//std::cout << "迭代时体积分数为" << vof_current << ",半径为" << ini_r << std::endl;
	}
	std::cout << "最终体积分数为" << vof_current << ",半径为" << ini_r << std::endl;
	if (ini_r < 0.8)
		return false;
	return true;
}

void Truss::genNormalTrussVoxelETH()
{
	double cell_length_half = cell_length / 2;	//模型大小的一半,与一些数值参数有关

	//重建范围
	Eigen::RowVector3d VminCut = { -1 * cell_length_half, -1 * cell_length_half, -1 * cell_length_half };
	Eigen::RowVector3d VmaxCut = { cell_length_half, cell_length_half, cell_length_half };

	Eigen::RowVector3i res = (resolution * ((VmaxCut - VminCut) / (VmaxCut - VminCut).maxCoeff())).cast<int>();

	grid.resize(res(0) * res(1) * res(2), 3);

	//构建采样点的坐标
	const auto lerpCut = [&](const int di, const int d)->double
	{return VminCut(d) + (double)di / (double)(res(d) - 1) * (VmaxCut(d) - VminCut(d)); };

#pragma omp parallel for
	for (int zi = 0; zi < res(2); zi++)
	{
		const double zc = lerpCut(zi, 2);
		for (int yi = 0; yi < res(1); yi++)
		{
			const double yc = lerpCut(yi, 1);
			for (int xi = 0; xi < res(0); xi++)
			{
				const double xc = lerpCut(xi, 0);
				grid.row(xi + res(0) * (yi + res(1) * zi)) = Eigen::RowVector3d(xc, yc, zc);
			}
		}
	}

	std::vector<Edge> edge_list = base_edges;
	//std::vector<Edge> edge_list = edges;


	//计算距离场的公式
	//看PPL的计算
	//将采样点xyz映射到子单元1/8或者1/48中,然后计算距离,只需要关照子单元中的几条杆
	const auto& function = [edge_list, cell_length_half](double x, double y, double z)->double {
		//超出边界,赋一个比较大的正值
		if (abs(x) > cell_length_half || abs(y) > cell_length_half || abs(z) > cell_length_half)
		{
			return 100;
		}
		//将采样点映射到子单元
		Eigen::Vector3d p = { x,y,z };
		getETHSymmetryPoint(p);

		double min_dis = 100;
		for (int i = 0; i < edge_list.size(); i++)
		{
			Eigen::Vector3d v1 = edge_list[i].v.row(0), v2 = edge_list[i].v.row(1);
			double dis = pointToSegDist(v1, v2, p) - edge_list[i].r_outer;
			if (dis < min_dis)
				min_dis = dis;
		}
		return min_dis;
	};
	voxel.resize(grid.rows());
	voxel.setZero();
#pragma omp parallel for
	for (int i = 0; i < grid.rows(); i++) {
		double indicator = function(grid(i, 0), grid(i, 1), grid(i, 2));
		if (indicator <= 0)
		{
			voxel(i) = 1;
		}
		else
			voxel(i) = 0;
	}
}

void Truss::genNormalTrussMeshETH()
{
	double cell_length_half = cell_length / 2;	//模型大小的一半,与一些数值参数有关

	//重建范围
	Eigen::RowVector3d Vmin = { -1 * cell_length_half - 0.1, -1 * cell_length_half - 0.1, -1 * cell_length_half - 0.1 };
	Eigen::RowVector3d Vmax = { cell_length_half + 0.1, cell_length_half + 0.1, cell_length_half + 0.1 };

	Eigen::RowVector3i res_mesh = (resolution_mesh * ((Vmax - Vmin) / (Vmax - Vmin).maxCoeff())).cast<int>();

	Eigen::MatrixXd GV(res_mesh(0) * res_mesh(1) * res_mesh(2), 3);

	//构建采样点的坐标
	const auto lerp = [&](const int di, const int d)->double
	{return Vmin(d) + (double)di / (double)(res_mesh(d) - 1) * (Vmax(d) - Vmin(d)); };

#pragma omp parallel for
	for (int zi = 0; zi < res_mesh(2); zi++)
	{
		const double z = lerp(zi, 2);
		for (int yi = 0; yi < res_mesh(1); yi++)
		{
			const double y = lerp(yi, 1);
			for (int xi = 0; xi < res_mesh(0); xi++)
			{
				const double x = lerp(xi, 0);
				GV.row(xi + res_mesh(0) * (yi + res_mesh(1) * zi)) = Eigen::RowVector3d(x, y, z);
			}
		}
	}

	std::vector<Edge> edge_list = base_edges;
	//std::vector<Edge> edge_list = edges;

	//计算距离场的公式
	//看PPL的计算
	//将采样点xyz映射到子单元1/8或者1/48中,然后计算距离,只需要关照子单元中的几条杆
	const auto& function = [edge_list, cell_length_half](double x, double y, double z)->double {
		//超出边界,赋一个比较大的正值
		if (abs(x) > cell_length_half || abs(y) > cell_length_half || abs(z) > cell_length_half)
		{
			return 100;
		}

		//将采样点映射到子单元
		Eigen::Vector3d p = { x,y,z };
		getETHSymmetryPoint(p);

		double min_dis = 100;
		for (int i = 0; i < edge_list.size(); i++)
		{
			Eigen::Vector3d v1 = edge_list[i].v.row(0), v2 = edge_list[i].v.row(1);
			double dis = pointToSegDist(v1, v2, p) - edge_list[i].r_outer;
			if (dis < min_dis)
				min_dis = dis;
		}
		return min_dis;
	};

	Eigen::VectorXd B(GV.rows());
#pragma omp parallel for
	for (int i = 0; i < GV.rows(); i++) {
		double indicator = function(GV(i, 0), GV(i, 1), GV(i, 2));
		B(i) = indicator;
	}
	igl::copyleft::marching_cubes(B, GV, res_mesh(0), res_mesh(1), res_mesh(2), V, F);

	//平移中心到原点
	if (V.rows() == 0 || !V.allFinite())
		throw std::runtime_error("No finite isosurface; check the full-cell skeleton, radius and mesh resolution");
	Eigen::Vector3d minV = V.colwise().minCoeff();
	Eigen::Vector3d maxV = V.colwise().maxCoeff();

	Eigen::Vector3d move = minV + (maxV - minV) / 2.0;
	for (int j = 0; j < V.rows(); j++)
	{
		V(j, 0) -= move(0);
		V(j, 1) -= move(1);
		V(j, 2) -= move(2);
	}

	//缩放为[-1,-1,-1]到[1,1,1]的边长为2的立方体
	Eigen::Vector3d factor = (maxV - minV) / 2.0;
	double xfac = (1.0 / factor(0));
	double yfac = (1.0 / factor(1));
	double zfac = (1.0 / factor(2));
	for (int j = 0; j < V.rows(); j++)
	{

		V(j, 0) *= xfac;
		V(j, 1) *= yfac;
		V(j, 2) *= zfac;
	}
}



/*void Truss::setHollowVolumeFraction(double vof)
{
	//计算当前体积分数
	genHollowTrussVoxel();
	double vof_1 = (double)voxel.sum() / (double)voxel.rows();
	//std::cout << "初始体积分数为" << vof_1 << ",半径为" << r_outer << std::endl;

	//将d_outer加一点点,然后计算斜率
	setROuter(r_outer + 0.5);
	genHollowTrussVoxel();
	double vof_2 = (double)voxel.sum() / (double)voxel.rows();
	//std::cout << "初始体积分数为" << vof_2 << ",半径为" << r_outer << std::endl;

	double k = ((vof_2 - vof_1) * 1.0) / 0.5;//斜率
	if (!std::isfinite(k)) return false;
	double ini_r = std::abs(k) < 1e-14 ? r_outer : ((vof - vof_2) / k + r_outer);//粗略计算目标相对密度的半径

	//迭代求解精确的半径,确保体积分数在+-0.5%
	setROuter(ini_r);
	genHollowTrussVoxel();
	double vof_current = (double)voxel.sum() / (double)voxel.rows();
	//std::cout << "迭代开始时体积分数为" << vof_current << ",半径为" << ini_r << std::endl;

	double factor = 1.0;	//步进长度因子
	int fit_iterations = 0;
	while (std::abs(vof - vof_current) > 0.005)
	{
		if (factor < 1e-10 || ++fit_iterations > 1000) return false;
		if (vof_current < vof)
		{
			//一点点增加
			while (vof - vof_current > 0.005)
			{
				if (++fit_iterations > 1000) return false;
				ini_r += 0.1 * factor;
				setROuter(ini_r);
				genHollowTrussVoxel();
				vof_current = (double)voxel.sum() / (double)voxel.rows();
			}
		}
		else if (vof_current > vof)
		{
			//一点点减少
			while (vof - vof_current < -0.005)
			{
				if (++fit_iterations > 1000) return false;
				ini_r -= 0.1 * factor;
				setROuter(ini_r);
				genHollowTrussVoxel();
				vof_current = (double)voxel.sum() / (double)voxel.rows();
			}
		}
		else {}
		factor = factor * 0.1;
		//std::cout << "迭代时体积分数为" << vof_current << ",半径为" << ini_r << std::endl;
	}
	std::cout << "最终体积分数为" << vof_current << ",半径为" << ini_r << std::endl;
}*/

/*void Truss::genHollowTrussVoxel()
{
	double cell_length_half = cell_length / 2;	//模型大小的一半,与一些数值参数有关

	//重建范围
	Eigen::RowVector3d VminCut = { -1 * cell_length_half, -1 * cell_length_half, -1 * cell_length_half };
	Eigen::RowVector3d VmaxCut = { cell_length_half, cell_length_half, cell_length_half };

	Eigen::RowVector3i res = (resolution * ((VmaxCut - VminCut) / (VmaxCut - VminCut).maxCoeff())).cast<int>();

	grid.resize(res(0) * res(1) * res(2), 3);

	//构建采样点的坐标
	const auto lerpCut = [&](const int di, const int d)->double
	{return VminCut(d) + (double)di / (double)(res(d) - 1) * (VmaxCut(d) - VminCut(d)); };

#pragma omp parallel for
	for (int zi = 0; zi < res(2); zi++)
	{
		const double zc = lerpCut(zi, 2);
		for (int yi = 0; yi < res(1); yi++)
		{
			const double yc = lerpCut(yi, 1);
			for (int xi = 0; xi < res(0); xi++)
			{
				const double xc = lerpCut(xi, 0);
				grid.row(xi + res(0) * (yi + res(1) * zi)) = Eigen::RowVector3d(xc, yc, zc);
			}
		}
	}

	//std::vector<Edge> edge_list = base_edges;
	std::vector<Edge> edge_list = edges;
	std::vector<Joint> joint_list = base_joints;
	double smf = smooth_factor;

	//计算距离场的公式
	//看PPL的计算
	//将采样点xyz映射到子单元1/8或者1/48中,然后计算距离,只需要关照子单元中的几条杆
	//返回值为pair<outer,inner>,outer<=0时体素为1;inner取反,所以inner也是<=0时体素为1
	const auto& function = [joint_list, edge_list, cell_length_half, smf](double x, double y, double z)->std::pair<double,double> {
		//超出边界,赋一个比较大的正值
		if (abs(x) > cell_length_half || abs(y) > cell_length_half || abs(z) > cell_length_half)
		{
			return std::make_pair(100, 100);
		}
		//将采样点映射到子单元
		// //初始1/48四面体的范围:
		//V1: p1=(0,0,0) p2=(1,0,0) p3=(1,0,1) p4=(1,1,1) 特征:x>=y, x>=z, y<=z
		//V2	: p1=(0,0,0) p2=(0,0,1) p3=(1,0,1) p4=(1,1,1) 特征:x>=y, x<=z, y<=z
		//V3: p1=(0,0,0) p2=(0,0,1) p3=(0,1,1) p4=(1,1,1) 特征:x<=y, x<=z, y<=z		 
		//其他四面体首先对称到xyz都为正数,然后如果y>z就按y=z平面对称上来,此时y<=z
		//这时候在三个四面体内,V1和V2关于x=z对称,V2和V3关于x=y对称
		//分别识别是否在V3,是否在V2,从而最后对称到V1中
		//todo:如果是子单元1/8,可能需要修改
		x = abs(x);
		y = abs(y);
		z = abs(z);
		if (y > z)
		{
			double temp = y;
			y = z;
			z = temp;
		}
		if (x < y)
		{
			double temp = x;
			x = y;
			y = temp;
		}
		if (x < z)
		{
			double temp = x;
			x = z;
			z = temp;
		}

		Eigen::Vector3d p = { x,y,z };
		double result_outer = 0, result_inner = 0;		//KS函数
		for (int i = 0; i < edge_list.size(); i++)
		{
			Eigen::Vector3d v1 = edge_list[i].v.row(0), v2 = edge_list[i].v.row(1);
			double dis_outer = pointToSegDist(v1, v2, p) - edge_list[i].r_outer;
			result_outer += exp(-1 / smf * dis_outer);	//smf是公式中的1/s_j
			double dis_inner = pointToSegDist(v1, v2, p) - edge_list[i].r_inner;
			result_inner += exp(-1 / smf * dis_inner);
		}
		result_outer = -smf * log(result_outer);
		result_inner = -smf * log(result_inner);
		return std::make_pair(result_outer, result_inner);
	};
	voxel.resize(grid.rows());
	voxel.setZero();
	for (int i = 0; i < grid.rows(); i++) {
		auto indicator = function(grid(i, 0), grid(i, 1), grid(i, 2));
		double outer = indicator.first, inner = -1.0 * indicator.second;		//注意inner取反了
		if (outer <= 0 && inner <= 0)
		{
			voxel(i) = 1;
		}
		else
			voxel(i) = 0;
	}
}*/

/*void Truss::genHollowTrussMesh()
{
	double cell_length_half = cell_length / 2;	//模型大小的一半,与一些数值参数有关

	//重建范围
	Eigen::RowVector3d Vmin = { -1 * cell_length_half - 0.1, -1 * cell_length_half - 0.1, -1 * cell_length_half - 0.1 };
	Eigen::RowVector3d Vmax = { cell_length_half + 0.1, cell_length_half + 0.1, cell_length_half + 0.1 };

	Eigen::RowVector3i res_mesh = (resolution_mesh * ((Vmax - Vmin) / (Vmax - Vmin).maxCoeff())).cast<int>();

	Eigen::MatrixXd GV(res_mesh(0) * res_mesh(1) * res_mesh(2), 3);

	//构建采样点的坐标
	const auto lerp = [&](const int di, const int d)->double
	{return Vmin(d) + (double)di / (double)(res_mesh(d) - 1) * (Vmax(d) - Vmin(d)); };

#pragma omp parallel for
	for (int zi = 0; zi < res_mesh(2); zi++)
	{
		const double z = lerp(zi, 2);
		for (int yi = 0; yi < res_mesh(1); yi++)
		{
			const double y = lerp(yi, 1);
			for (int xi = 0; xi < res_mesh(0); xi++)
			{
				const double x = lerp(xi, 0);
				GV.row(xi + res_mesh(0) * (yi + res_mesh(1) * zi)) = Eigen::RowVector3d(x, y, z);
			}
		}
	}

	//std::vector<Edge> edge_list = base_edges;
	std::vector<Edge> edge_list = edges;
	std::vector<Joint> joint_list = base_joints;
	double smf = smooth_factor;

	//计算距离场的公式
	//看PPL的计算
	//将采样点xyz映射到子单元1/8或者1/48中,然后计算距离,只需要关照子单元中的几条杆
	const auto& function = [joint_list, edge_list, cell_length_half, smf](double x, double y, double z)->std::pair<double,double> {
		//超出边界,赋一个比较大的正值
		if (abs(x) > cell_length_half || abs(y) > cell_length_half || abs(z) > cell_length_half)
		{
			return std::make_pair(100, 100);
		}
		//将采样点映射到子单元
		// //初始1/48四面体的范围:
		//V1: p1=(0,0,0) p2=(1,0,0) p3=(1,0,1) p4=(1,1,1) 特征:x>=y, x>=z, y<=z
		//V2	: p1=(0,0,0) p2=(0,0,1) p3=(1,0,1) p4=(1,1,1) 特征:x>=y, x<=z, y<=z
		//V3: p1=(0,0,0) p2=(0,0,1) p3=(0,1,1) p4=(1,1,1) 特征:x<=y, x<=z, y<=z		 
		//其他四面体首先对称到xyz都为正数,然后如果y>z就按y=z平面对称上来,此时y<=z
		//这时候在三个四面体内,V1和V2关于x=z对称,V2和V3关于x=y对称
		//分别识别是否在V3,是否在V2,从而最后对称到V1中
		//todo:如果是子单元1/8,可能需要修改
		x = abs(x);
		y = abs(y);
		z = abs(z);
		if (y > z)
		{
			double temp = y;
			y = z;
			z = temp;
		}
		if (x < y)
		{
			double temp = x;
			x = y;
			y = temp;
		}
		if (x < z)
		{
			double temp = x;
			x = z;
			z = temp;
		}

		Eigen::Vector3d p = { x,y,z };
		double result_outer = 0, result_inner = 0;		//KS函数
		for (int i = 0; i < edge_list.size(); i++)
		{
			Eigen::Vector3d v1 = edge_list[i].v.row(0), v2 = edge_list[i].v.row(1);
			double dis_outer = pointToSegDist(v1, v2, p) - edge_list[i].r_outer;
			result_outer += exp(-1 / smf * dis_outer);	//smf是公式中的1/s_j
			double dis_inner = pointToSegDist(v1, v2, p) - edge_list[i].r_inner;
			result_inner += exp(-1 / smf * dis_inner);
		}
		result_outer = -smf * log(result_outer);
		result_inner = -smf * log(result_inner);
		return std::make_pair(result_outer, result_inner);
	};

	Eigen::VectorXd B(GV.rows());
#pragma omp parallel for
	for (int i = 0; i < GV.rows(); i++) {
		auto indicator = function(GV(i, 0), GV(i, 1), GV(i, 2));
		double outer = indicator.first, inner = -1.0 * indicator.second;		//注意inner取反了
		B(i) = std::max(outer, inner);
	}
	igl::copyleft::marching_cubes(B, GV, res_mesh(0), res_mesh(1), res_mesh(2), V, F);

	//平移中心到原点
	if (V.rows() == 0 || !V.allFinite())
		throw std::runtime_error("No finite isosurface; check the full-cell skeleton, radius and mesh resolution");
	Eigen::Vector3d minV = V.colwise().minCoeff();
	Eigen::Vector3d maxV = V.colwise().maxCoeff();

	Eigen::Vector3d move = minV + (maxV - minV) / 2.0;
	for (int j = 0; j < V.rows(); j++)
	{
		V(j, 0) -= move(0);
		V(j, 1) -= move(1);
		V(j, 2) -= move(2);
	}

	//缩放为[-1,-1,-1]到[1,1,1]的边长为2的立方体
	Eigen::Vector3d factor = (maxV - minV) / 2.0;
	double xfac = (1.0 / factor(0));
	double yfac = (1.0 / factor(1));
	double zfac = (1.0 / factor(2));
	for (int j = 0; j < V.rows(); j++)
	{

		V(j, 0) *= xfac;
		V(j, 1) *= yfac;
		V(j, 2) *= zfac;
	}
}*/


bool Truss::isInnerVoxel(int index)
{
	int range = resolution;
	int z = index / (range * range);
	int y = index % (range * range) / range;
	int x = index % (range * range) % range;
	int direction[6] = { index + 1,index - 1,index + range,index - range,index + range * range,index - range * range };
	Eigen::Vector3i coor = { x,y,z };

	//在范围内且邻居体素为空,则不是内部体素
	if (x > 0 && x < range - 1 && y>0 && y < range - 1 && z>0 && z < range - 1) {
		for (int i = 0; i < 6; i++) {
			if (voxel[direction[i]] == 0) {
				return false;
			}
		}
	}
	else {	//处理单胞边界上的体素
		for (int i = 0; i < 3; i++) {
			Eigen::Vector3i neibor1 = coor;
			Eigen::Vector3i neibor2 = coor;
			if (coor(i) == 0) {
				neibor1(i) += 1;
				if (voxel[neibor1(0) + range * (neibor1(1) + range * (neibor1(2)))] == 0) return false;
			}
			else if (coor(i) == range - 1) {
				neibor2(i) -= 1;
				if (voxel[neibor2(0) + range * (neibor2(1) + range * (neibor2(2)))] == 0) return false;
			}
			else {
				neibor1(i)++;
				neibor2(i)--;
				if (voxel[neibor1(0) + range * (neibor1(1) + range * (neibor1(2)))] == 0) return false;
				if (voxel[neibor2(0) + range * (neibor2(1) + range * (neibor2(2)))] == 0) return false;
			}
		}
	}
	return true;
}


void Truss::getTargetCoors(Eigen::MatrixXd& tgrid, Eigen::VectorXi& mapping)
{
	double cell_length_half = cell_length / 2;	//模型大小的一半,与一些数值参数有关

	//重建范围
	Eigen::RowVector3d VminCut = { -1 * cell_length_half, -1 * cell_length_half, -1 * cell_length_half };
	Eigen::RowVector3d VmaxCut = { cell_length_half, cell_length_half, cell_length_half };

	Eigen::RowVector3i res = (resolution * ((VmaxCut - VminCut) / (VmaxCut - VminCut).maxCoeff())).cast<int>();

	grid.resize(res(0) * res(1) * res(2), 3);

	//构建采样点的坐标
	const auto lerpCut = [&](const int di, const int d)->double
	{return VminCut(d) + (double)di / (double)(res(d) - 1) * (VmaxCut(d) - VminCut(d)); };

#pragma omp parallel for
	for (int zi = 0; zi < res(2); zi++)
	{
		const double zc = lerpCut(zi, 2);
		for (int yi = 0; yi < res(1); yi++)
		{
			const double yc = lerpCut(yi, 1);
			for (int xi = 0; xi < res(0); xi++)
			{
				const double xc = lerpCut(xi, 0);
				grid.row(xi + res(0) * (yi + res(1) * zi)) = Eigen::RowVector3d(xc, yc, zc);
			}
		}
	}
	//取出值为1的体素坐标放到tgrid中
	int tgrid_index = 0;
	tgrid.resize(voxel.sum(), 3);
	mapping.resize(grid.rows());
	for (int i = 0; i < grid.rows(); i++)
	{
		if (voxel(i) == 1) {
			tgrid.row(tgrid_index) = grid.row(i);
			mapping(i) = tgrid_index;
			tgrid_index++;
		}
		else {
			mapping(i) = -1;
		}
	}

}


void Truss::genTrussGrid()
{
	double cell_length_half = cell_length / 2;	//模型大小的一半,与一些数值参数有关

	//重建范围
	Eigen::RowVector3d VminCut = { -1 * cell_length_half, -1 * cell_length_half, -1 * cell_length_half };
	Eigen::RowVector3d VmaxCut = { cell_length_half, cell_length_half, cell_length_half };

	Eigen::RowVector3i res = (resolution * ((VmaxCut - VminCut) / (VmaxCut - VminCut).maxCoeff())).cast<int>();

	grid.resize((res(0) + 1) * (res(1) + 1) * (res(2) + 1), 3);
	grid.setZero();

	//构建采样点的坐标
	const auto lerpCut = [&](const int di, const int d)->double
	{return VminCut(d) + (double)di / (double)(res(d) - 1) * (VmaxCut(d) - VminCut(d)); };

#pragma omp parallel for
	for (int zi = 0; zi <= res(2); zi++)
	{
		const double zc = lerpCut(zi, 2);
		for (int yi = 0; yi <= res(1); yi++)
		{
			const double yc = lerpCut(yi, 1);
			for (int xi = 0; xi <= res(0); xi++)
			{
				const double xc = lerpCut(xi, 0);
				grid.row(xi + (res(0) + 1) * (yi + (res(1) + 1) * zi)) = Eigen::RowVector3d(xc, yc, zc);
			}
		}
	}
}



void Truss::genSolidTargetSdf(Eigen::MatrixXd& tgrid)
{
	double cell_length_half = cell_length / 2;	//模型大小的一半,与一些数值参数有关

	//重建范围
	Eigen::RowVector3d VminCut = { -1 * cell_length_half, -1 * cell_length_half, -1 * cell_length_half };
	Eigen::RowVector3d VmaxCut = { cell_length_half, cell_length_half, cell_length_half };

	Eigen::RowVector3i res = (resolution * ((VmaxCut - VminCut) / (VmaxCut - VminCut).maxCoeff())).cast<int>();

	grid.resize(res(0) * res(1) * res(2), 3);

	//构建采样点的坐标
	const auto lerpCut = [&](const int di, const int d)->double
	{return VminCut(d) + (double)di / (double)(res(d) - 1) * (VmaxCut(d) - VminCut(d)); };

#pragma omp parallel for
	for (int zi = 0; zi < res(2); zi++)
	{
		const double zc = lerpCut(zi, 2);
		for (int yi = 0; yi < res(1); yi++)
		{
			const double yc = lerpCut(yi, 1);
			for (int xi = 0; xi < res(0); xi++)
			{
				const double xc = lerpCut(xi, 0);
				grid.row(xi + res(0) * (yi + res(1) * zi)) = Eigen::RowVector3d(xc, yc, zc);
			}
		}
	}

	//std::vector<Edge> edge_list = base_edges;
	std::vector<Edge> edge_list = edges;
	std::vector<Joint> joint_list = nodes;
	//double smf = smooth_factor;

	//计算距离场的公式
	//看PPL的计算
	//将采样点xyz映射到子单元1/8或者1/48中,然后计算距离,只需要关照子单元中的几条杆
	const auto& function = [edge_list, joint_list, cell_length_half](double x, double y, double z)->double {
		//超出边界,赋一个比较大的正值
		if (abs(x) > cell_length_half || abs(y) > cell_length_half || abs(z) > cell_length_half)
		{
			return 100;
		}
		//将采样点映射到子单元
		Eigen::Vector3d p = { x,y,z };
		getPanettaSymmetryPoint(p);

		std::vector<double> result;	//存放每个joint的KS函数值,到最后再套一层总的KS函数

		for (int i = 0; i < joint_list.size(); i++)
		{
			double result_joint = 0;		//当前joint的KS函数值
			double smf = joint_list[i].smooth_factor;
			for (int j = 0; j < joint_list[i].associated_edges.size(); j++)
			{
				int index = joint_list[i].associated_edges[j];	//相邻边在edge_list中的索引
				Eigen::Vector3d v1 = edge_list[index].v.row(0), v2 = edge_list[index].v.row(1);
				double dis = pointToSegDist(v1, v2, p) - edge_list[index].r_outer;
				result_joint += exp(-1 / smf * dis);	//smf是公式中的1/s_j
			}
			result_joint = -smf * log(result_joint);
			result.push_back(result_joint);
		}

		double sum_KS = 0;	//总的KS函数值
		double sum_smf = 0.5;
		for (int i = 0; i < result.size(); i++)
		{
			sum_KS += exp(-1 / sum_smf * result[i]);
		}
		sum_KS = -sum_smf * log(sum_KS);
		return sum_KS;
	};

	outer_sdf.setZero(tgrid.rows());
	//inner_sdf.setZero(tgrid.rows());

#pragma omp parallel for
	for (int i = 0; i < tgrid.rows(); i++)
	{
		double indicator = function(tgrid(i, 0), tgrid(i, 1), tgrid(i, 2));
		outer_sdf(i) = indicator;
		//实心杆不用inner_sdf
	}

	voxel.resize(grid.rows());
	voxel.setZero();
#pragma omp parallel for
	for (int i = 0; i < grid.rows(); i++) {
		double indicator = function(grid(i, 0), grid(i, 1), grid(i, 2));
		if (indicator <= 0)
		{
			voxel(i) = 1;
		}
		else
			voxel(i) = 0;
	}

	tvol = voxel.sum() / (double)voxel.rows();
}


void Truss::getSolidVelocities(std::vector<Eigen::VectorXd>& velosities)
{
	//目前设置delta=0.7 可以自行调整
	double delta = 0.7 * (grid.row(1) - grid.row(0)).x();
	double delta2 = pow(delta, 2);	//没用上
	delta = delta2;
	//double coef = 0.75 / delta;		//没用上

	Truss PTL2 = *this;
	double perturbation_R = 1.0, perturbation_smf = 0.3;		//步进值,需要自行调整

	//计算当前的体素和距离场
	Eigen::MatrixXd tgrid;	//值为1的体素坐标
	Eigen::VectorXi map;
	getTargetCoors(tgrid, map);
	genSolidTargetSdf(tgrid);
	int vsize = tgrid.rows();

	velosities.clear();
	dV.clear();
	dH.clear();

	//计算R改变后的体素和距离场,并计算差值
	for (int index = 0; index < base_edges.size(); index++)
	{
		//int index = 0;
		Eigen::VectorXd vv;
		double current_r_outer = PTL2.getROuter(index);
		PTL2.setROuter(current_r_outer + perturbation_R, index);
		PTL2.genSolidTargetSdf(tgrid);
		PTL2.setROuter(current_r_outer, index);	//复原

		vv.resize(vsize);
		vv.setZero();
		Eigen::VectorXd sdf_diff;	//两个距离场的差
		sdf_diff.setZero(vsize);
		Eigen::VectorXd outer_sdf_diff = outer_sdf - PTL2.outer_sdf;
		//Eigen::VectorXd inner_sdf_diff = inner_sdf - s2.inner_sdf;

		Eigen::VectorXi vdiff = PTL2.voxel - voxel;	//为什么这里是PTL2-PTL1? 因为距离场diff取了绝对值

		Eigen::VectorXd diffp(vsize);
		Eigen::VectorXd diffn(vsize);
		Eigen::VectorXi vdiffp(voxel.rows());
		Eigen::VectorXi vdiffn(voxel.rows());
		diffp.setZero();
		diffn.setZero();
		vdiffp.setZero();
		vdiffn.setZero();

#pragma omp parallel for
		for (int j = 0; j < voxel.rows(); j++) {
			if (voxel(j) == 1) {
				//需要是边界上的体素点
				if ((!isInnerVoxel(j))) {
					int k = map(j);	//取出在tgrid中的坐标

					//Dirac function
					//identify outer boundary
					if (abs(outer_sdf(k)) <= delta) {
						sdf_diff(k) = outer_sdf_diff(k);
					}
					//inner boundary
					//else if (abs(inner_sdf(k) <= delta)) {
					//	sdf_diff(k) = inner_sdf_diff(k);
					//}
					else {
						sdf_diff(k) = 0;
					}
					if (sdf_diff(k) > 0) {
						diffp(k) = sdf_diff(k);
					}
					if (sdf_diff(k) < 0) {
						diffn(k) = sdf_diff(k);
					}
				}
			}
			if (vdiff(j) > 0) {
				vdiffp(j) = vdiff(j);
			}
			else if (vdiff(j) < 0) {
				vdiffn(j) = vdiff(j);
			}
		}

		double sum_diffp = diffp.sum();
		double sum_diffn = diffn.sum();
		double sum_vdiffp = vdiffp.sum();
		double sum_vdiffn = vdiffn.sum();

		double factor1, factor2;
		if (sum_diffp == 0) {
			factor1 = 0;
		}
		else {
			factor1 = sum_vdiffp / sum_diffp;
		}
		if (sum_diffn == 0) {

			factor2 = 0;
		}
		else {
			factor2 = sum_vdiffn / sum_diffn;
		}

		vv = factor1 * diffp + factor2 * diffn;
		vv = vv / perturbation_R;

		double dv = (PTL2.tvol - tvol) / perturbation_R;
		velosities.push_back(vv);
		dV.push_back(dv);
	}
	//计算smf改变后的体素和距离场,并计算差值
	for (int index = 0; index < base_joints.size(); index++)
	{
		Eigen::VectorXd vv;
		double current_smooth_factor = PTL2.getSmoothFactor(index);
		PTL2.setSmoothFactor(current_smooth_factor + perturbation_smf, index);
		PTL2.genSolidTargetSdf(tgrid);
		PTL2.setSmoothFactor(current_smooth_factor, index);	//复原

		vv.resize(vsize);
		vv.setZero();
		Eigen::VectorXd sdf_diff;	//两个距离场的差
		sdf_diff.setZero(vsize);
		Eigen::VectorXd outer_sdf_diff = outer_sdf - PTL2.outer_sdf;
		//Eigen::VectorXd inner_sdf_diff = inner_sdf - s2.inner_sdf;

		Eigen::VectorXi vdiff = PTL2.voxel - voxel;	//为什么这里是PTL2-PTL1? 因为距离场diff取了绝对值

		Eigen::VectorXd diffp(vsize);
		Eigen::VectorXd diffn(vsize);
		Eigen::VectorXi vdiffp(voxel.rows());
		Eigen::VectorXi vdiffn(voxel.rows());
		diffp.setZero();
		diffn.setZero();
		vdiffp.setZero();
		vdiffn.setZero();

#pragma omp parallel for
		for (int j = 0; j < voxel.rows(); j++) {
			if (voxel(j) == 1) {
				//需要是边界上的体素点
				if ((!isInnerVoxel(j))) {
					int k = map(j);	//取出在tgrid中的坐标

					//Dirac function
					//identify outer boundary
					if (abs(outer_sdf(k)) <= delta) {
						sdf_diff(k) = outer_sdf_diff(k);
					}
					//inner boundary
					//else if (abs(inner_sdf(k) <= delta)) {
					//	sdf_diff(k) = inner_sdf_diff(k);
					//}
					else {
						sdf_diff(k) = 0;
					}
					if (sdf_diff(k) > 0) {
						diffp(k) = sdf_diff(k);
					}
					if (sdf_diff(k) < 0) {
						diffn(k) = sdf_diff(k);
					}
				}
			}
			if (vdiff(j) > 0) {
				vdiffp(j) = vdiff(j);
			}
			else if (vdiff(j) < 0) {
				vdiffn(j) = vdiff(j);
			}
		}

		double sum_diffp = diffp.sum();
		double sum_diffn = diffn.sum();
		double sum_vdiffp = vdiffp.sum();
		double sum_vdiffn = vdiffn.sum();

		double factor1, factor2;
		if (sum_diffp == 0) {
			factor1 = 0;
		}
		else {
			factor1 = sum_vdiffp / sum_diffp;
		}
		if (sum_diffn == 0) {

			factor2 = 0;
		}
		else {
			factor2 = sum_vdiffn / sum_diffn;
		}

		vv = factor1 * diffp + factor2 * diffn;
		vv = vv / perturbation_smf;

		double dv = (PTL2.tvol - tvol) / perturbation_smf;
		velosities.push_back(vv);
		dV.push_back(dv);
	}
}



void Truss::genSolidTargetSdfETH(Eigen::MatrixXd& tgrid)
{
	double cell_length_half = cell_length / 2;	//模型大小的一半,与一些数值参数有关

	//重建范围
	Eigen::RowVector3d VminCut = { -1 * cell_length_half, -1 * cell_length_half, -1 * cell_length_half };
	Eigen::RowVector3d VmaxCut = { cell_length_half, cell_length_half, cell_length_half };

	Eigen::RowVector3i res = (resolution * ((VmaxCut - VminCut) / (VmaxCut - VminCut).maxCoeff())).cast<int>();

	grid.resize(res(0) * res(1) * res(2), 3);

	//构建采样点的坐标
	const auto lerpCut = [&](const int di, const int d)->double
	{return VminCut(d) + (double)di / (double)(res(d) - 1) * (VmaxCut(d) - VminCut(d)); };

#pragma omp parallel for
	for (int zi = 0; zi < res(2); zi++)
	{
		const double zc = lerpCut(zi, 2);
		for (int yi = 0; yi < res(1); yi++)
		{
			const double yc = lerpCut(yi, 1);
			for (int xi = 0; xi < res(0); xi++)
			{
				const double xc = lerpCut(xi, 0);
				grid.row(xi + res(0) * (yi + res(1) * zi)) = Eigen::RowVector3d(xc, yc, zc);
			}
		}
	}

	//std::vector<Edge> edge_list = base_edges;
	std::vector<Edge> edge_list = edges;
	std::vector<Joint> joint_list = nodes;
	//double smf = smooth_factor;

	//计算距离场的公式
	//看PPL的计算
	//将采样点xyz映射到子单元1/8或者1/48中,然后计算距离,只需要关照子单元中的几条杆
	const auto& function = [edge_list, joint_list, cell_length_half](double x, double y, double z)->double {
		//超出边界,赋一个比较大的正值
		if (abs(x) > cell_length_half || abs(y) > cell_length_half || abs(z) > cell_length_half)
		{
			return 100;
		}
		//将采样点映射到子单元
		Eigen::Vector3d p = { x,y,z };
		getETHSymmetryPoint(p);

		std::vector<double> result;	//存放每个joint的KS函数值,到最后再套一层总的KS函数

		for (int i = 0; i < joint_list.size(); i++)
		{
			double result_joint = 0;		//当前joint的KS函数值
			double smf = joint_list[i].smooth_factor;
			for (int j = 0; j < joint_list[i].associated_edges.size(); j++)
			{
				int index = joint_list[i].associated_edges[j];	//相邻边在edge_list中的索引
				Eigen::Vector3d v1 = edge_list[index].v.row(0), v2 = edge_list[index].v.row(1);
				double dis = pointToSegDist(v1, v2, p) - edge_list[index].r_outer;
				result_joint += exp(-1 / smf * dis);	//smf是公式中的1/s_j
			}
			result_joint = -smf * log(result_joint);
			result.push_back(result_joint);
		}

		double sum_KS = 0;	//总的KS函数值
		double sum_smf = 0.5;
		for (int i = 0; i < result.size(); i++)
		{
			sum_KS += exp(-1 / sum_smf * result[i]);
		}
		sum_KS = -sum_smf * log(sum_KS);
		return sum_KS;
	};

	outer_sdf.setZero(tgrid.rows());
	//inner_sdf.setZero(tgrid.rows());

#pragma omp parallel for
	for (int i = 0; i < tgrid.rows(); i++)
	{
		double indicator = function(tgrid(i, 0), tgrid(i, 1), tgrid(i, 2));
		outer_sdf(i) = indicator;
		//实心杆不用inner_sdf
	}

	voxel.resize(grid.rows());
	voxel.setZero();
#pragma omp parallel for
	for (int i = 0; i < grid.rows(); i++) {
		double indicator = function(grid(i, 0), grid(i, 1), grid(i, 2));
		if (indicator <= 0)
		{
			voxel(i) = 1;
		}
		else
			voxel(i) = 0;
	}

	tvol = voxel.sum() / (double)voxel.rows();
}


void Truss::getSolidVelocitiesETH(std::vector<Eigen::VectorXd>& velosities)
{
	//目前设置delta=0.7 可以自行调整
	double delta = 0.7 * (grid.row(1) - grid.row(0)).x();
	double delta2 = pow(delta, 2);	//没用上
	delta = delta2;
	//double coef = 0.75 / delta;		//没用上

	Truss PTL2 = *this;
	double perturbation_R = 1.0, perturbation_smf = 0.3;		//步进值,需要自行调整

	//计算当前的体素和距离场
	Eigen::MatrixXd tgrid;	//值为1的体素坐标
	Eigen::VectorXi map;
	getTargetCoors(tgrid, map);
	genSolidTargetSdfETH(tgrid);
	int vsize = tgrid.rows();

	velosities.clear();
	dV.clear();
	dH.clear();

	//计算R改变后的体素和距离场,并计算差值
	for (int index = 0; index < base_edges.size(); index++)
	{
		//int index = 0;
		Eigen::VectorXd vv;
		double current_r_outer = PTL2.getROuter(index);
		PTL2.setROuter(current_r_outer + perturbation_R, index);
		PTL2.genSolidTargetSdfETH(tgrid);
		PTL2.setROuter(current_r_outer, index);	//复原

		vv.resize(vsize);
		vv.setZero();
		Eigen::VectorXd sdf_diff;	//两个距离场的差
		sdf_diff.setZero(vsize);
		Eigen::VectorXd outer_sdf_diff = outer_sdf - PTL2.outer_sdf;
		//Eigen::VectorXd inner_sdf_diff = inner_sdf - s2.inner_sdf;

		Eigen::VectorXi vdiff = PTL2.voxel - voxel;	//为什么这里是PTL2-PTL1? 因为距离场diff取了绝对值

		Eigen::VectorXd diffp(vsize);
		Eigen::VectorXd diffn(vsize);
		Eigen::VectorXi vdiffp(voxel.rows());
		Eigen::VectorXi vdiffn(voxel.rows());
		diffp.setZero();
		diffn.setZero();
		vdiffp.setZero();
		vdiffn.setZero();

#pragma omp parallel for
		for (int j = 0; j < voxel.rows(); j++) {
			if (voxel(j) == 1) {
				//需要是边界上的体素点
				if ((!isInnerVoxel(j))) {
					int k = map(j);	//取出在tgrid中的坐标

					//Dirac function
					//identify outer boundary
					if (abs(outer_sdf(k)) <= delta) {
						sdf_diff(k) = outer_sdf_diff(k);
					}
					//inner boundary
					//else if (abs(inner_sdf(k) <= delta)) {
					//	sdf_diff(k) = inner_sdf_diff(k);
					//}
					else {
						sdf_diff(k) = 0;
					}
					if (sdf_diff(k) > 0) {
						diffp(k) = sdf_diff(k);
					}
					if (sdf_diff(k) < 0) {
						diffn(k) = sdf_diff(k);
					}
				}
			}
			if (vdiff(j) > 0) {
				vdiffp(j) = vdiff(j);
			}
			else if (vdiff(j) < 0) {
				vdiffn(j) = vdiff(j);
			}
		}

		double sum_diffp = diffp.sum();
		double sum_diffn = diffn.sum();
		double sum_vdiffp = vdiffp.sum();
		double sum_vdiffn = vdiffn.sum();

		double factor1, factor2;
		if (sum_diffp == 0) {
			factor1 = 0;
		}
		else {
			factor1 = sum_vdiffp / sum_diffp;
		}
		if (sum_diffn == 0) {

			factor2 = 0;
		}
		else {
			factor2 = sum_vdiffn / sum_diffn;
		}

		vv = factor1 * diffp + factor2 * diffn;
		vv = vv / perturbation_R;

		double dv = (PTL2.tvol - tvol) / perturbation_R;
		velosities.push_back(vv);
		dV.push_back(dv);
	}
	//计算smf改变后的体素和距离场,并计算差值
	for (int index = 0; index < base_joints.size(); index++)
	{
		Eigen::VectorXd vv;
		double current_smooth_factor = PTL2.getSmoothFactor(index);
		PTL2.setSmoothFactor(current_smooth_factor + perturbation_smf, index);
		PTL2.genSolidTargetSdfETH(tgrid);
		PTL2.setSmoothFactor(current_smooth_factor, index);	//复原

		vv.resize(vsize);
		vv.setZero();
		Eigen::VectorXd sdf_diff;	//两个距离场的差
		sdf_diff.setZero(vsize);
		Eigen::VectorXd outer_sdf_diff = outer_sdf - PTL2.outer_sdf;
		//Eigen::VectorXd inner_sdf_diff = inner_sdf - s2.inner_sdf;

		Eigen::VectorXi vdiff = PTL2.voxel - voxel;	//为什么这里是PTL2-PTL1? 因为距离场diff取了绝对值

		Eigen::VectorXd diffp(vsize);
		Eigen::VectorXd diffn(vsize);
		Eigen::VectorXi vdiffp(voxel.rows());
		Eigen::VectorXi vdiffn(voxel.rows());
		diffp.setZero();
		diffn.setZero();
		vdiffp.setZero();
		vdiffn.setZero();

#pragma omp parallel for
		for (int j = 0; j < voxel.rows(); j++) {
			if (voxel(j) == 1) {
				//需要是边界上的体素点
				if ((!isInnerVoxel(j))) {
					int k = map(j);	//取出在tgrid中的坐标

					//Dirac function
					//identify outer boundary
					if (abs(outer_sdf(k)) <= delta) {
						sdf_diff(k) = outer_sdf_diff(k);
					}
					//inner boundary
					//else if (abs(inner_sdf(k) <= delta)) {
					//	sdf_diff(k) = inner_sdf_diff(k);
					//}
					else {
						sdf_diff(k) = 0;
					}
					if (sdf_diff(k) > 0) {
						diffp(k) = sdf_diff(k);
					}
					if (sdf_diff(k) < 0) {
						diffn(k) = sdf_diff(k);
					}
				}
			}
			if (vdiff(j) > 0) {
				vdiffp(j) = vdiff(j);
			}
			else if (vdiff(j) < 0) {
				vdiffn(j) = vdiff(j);
			}
		}

		double sum_diffp = diffp.sum();
		double sum_diffn = diffn.sum();
		double sum_vdiffp = vdiffp.sum();
		double sum_vdiffn = vdiffn.sum();

		double factor1, factor2;
		if (sum_diffp == 0) {
			factor1 = 0;
		}
		else {
			factor1 = sum_vdiffp / sum_diffp;
		}
		if (sum_diffn == 0) {

			factor2 = 0;
		}
		else {
			factor2 = sum_vdiffn / sum_diffn;
		}

		vv = factor1 * diffp + factor2 * diffn;
		vv = vv / perturbation_smf;

		double dv = (PTL2.tvol - tvol) / perturbation_smf;
		velosities.push_back(vv);
		dV.push_back(dv);
	}
}


/*void Truss::genHollowTargetSdf(Eigen::MatrixXd& tgrid)
{
	double cell_length_half = cell_length / 2;	//模型大小的一半,与一些数值参数有关

	//重建范围
	Eigen::RowVector3d VminCut = { -1 * cell_length_half, -1 * cell_length_half, -1 * cell_length_half };
	Eigen::RowVector3d VmaxCut = { cell_length_half, cell_length_half, cell_length_half };

	Eigen::RowVector3i res = (resolution * ((VmaxCut - VminCut) / (VmaxCut - VminCut).maxCoeff())).cast<int>();

	grid.resize(res(0) * res(1) * res(2), 3);

	//构建采样点的坐标
	const auto lerpCut = [&](const int di, const int d)->double
	{return VminCut(d) + (double)di / (double)(res(d) - 1) * (VmaxCut(d) - VminCut(d)); };

#pragma omp parallel for
	for (int zi = 0; zi < res(2); zi++)
	{
		const double zc = lerpCut(zi, 2);
		for (int yi = 0; yi < res(1); yi++)
		{
			const double yc = lerpCut(yi, 1);
			for (int xi = 0; xi < res(0); xi++)
			{
				const double xc = lerpCut(xi, 0);
				grid.row(xi + res(0) * (yi + res(1) * zi)) = Eigen::RowVector3d(xc, yc, zc);
			}
		}
	}

	//std::vector<Edge> edge_list = base_edges;
	std::vector<Edge> edge_list = edges;
	std::vector<Joint> joint_list = base_joints;
	double smf = smooth_factor;

	//计算距离场的公式
	//看PPL的计算
	//将采样点xyz映射到子单元1/8或者1/48中,然后计算距离,只需要关照子单元中的几条杆
	const auto& function = [joint_list, edge_list, cell_length_half, smf](double x, double y, double z)->double {
		//超出边界,赋一个比较大的正值
		if (abs(x) > cell_length_half || abs(y) > cell_length_half || abs(z) > cell_length_half)
		{
			return 100;
		}
		//将采样点映射到子单元
		// //初始1/48四面体的范围:
		//V1: p1=(0,0,0) p2=(1,0,0) p3=(1,0,1) p4=(1,1,1) 特征:x>=y, x>=z, y<=z
		//V2	: p1=(0,0,0) p2=(0,0,1) p3=(1,0,1) p4=(1,1,1) 特征:x>=y, x<=z, y<=z
		//V3: p1=(0,0,0) p2=(0,0,1) p3=(0,1,1) p4=(1,1,1) 特征:x<=y, x<=z, y<=z		 
		//其他四面体首先对称到xyz都为正数,然后如果y>z就按y=z平面对称上来,此时y<=z
		//这时候在三个四面体内,V1和V2关于x=z对称,V2和V3关于x=y对称
		//分别识别是否在V3,是否在V2,从而最后对称到V1中
		//todo:如果是子单元1/8,可能需要修改
		x = abs(x);
		y = abs(y);
		z = abs(z);
		if (y > z)
		{
			double temp = y;
			y = z;
			z = temp;
		}
		if (x < y)
		{
			double temp = x;
			x = y;
			y = temp;
		}
		if (x < z)
		{
			double temp = x;
			x = z;
			z = temp;
		}

		Eigen::Vector3d p = { x,y,z };
		double result = 0;	//KS函数
		for (int i = 0; i < edge_list.size(); i++)
		{
			Eigen::Vector3d v1 = edge_list[i].v.row(0), v2 = edge_list[i].v.row(1);
			double dis = pointToSegDist(v1, v2, p) - edge_list[i].r_outer;
			result += exp(-1 / smf * dis);	//smf是公式中的1/s_j
		}
		result = -smf * log(result);
		return result;
	};


	outer_sdf.setZero(tgrid.rows());
	//inner_sdf.setZero(tgrid.rows());

#pragma omp parallel for
	for (int i = 0; i < tgrid.rows(); i++)
	{
		double indicator = function(tgrid(i, 0), tgrid(i, 1), tgrid(i, 2));
		outer_sdf(i) = indicator;
		//实心杆不用inner_sdf
	}

	voxel.resize(grid.rows());
	voxel.setZero();
#pragma omp parallel for
	for (int i = 0; i < grid.rows(); i++) {
		double indicator = function(grid(i, 0), grid(i, 1), grid(i, 2));
		if (indicator <= 0)
		{
			voxel(i) = 1;
		}
		else
			voxel(i) = 0;
	}

	tvol = voxel.sum() / (double)voxel.rows();
}*/

/*void Truss::getHollowVelocities(std::vector<Eigen::VectorXd>& velosities)
{

}*/


bool Truss::setFusionVolumeFraction(double vof)
{
	//计算当前体积分数
	genFusionTrussVoxel();
	double vof_1 = (double)voxel.sum() / (double)voxel.rows();
	//std::cout << "初始体积分数为" << vof_1 << ",半径为" << r_outer << std::endl;

	//将d_outer加一点点,然后计算斜率
	setROuter(r_outer + 0.5, -1);
	genFusionTrussVoxel();
	double vof_2 = (double)voxel.sum() / (double)voxel.rows();
	//std::cout << "初始体积分数为" << vof_2 << ",半径为" << r_outer << std::endl;

	double k = ((vof_2 - vof_1) * 1.0) / 0.5;//斜率
	if (!std::isfinite(k)) return false;
	double ini_r = std::abs(k) < 1e-14 ? r_outer : ((vof - vof_2) / k + r_outer);//粗略计算目标相对密度的半径
	if (ini_r < 0)
		ini_r = 1.0;		//从最小开始

	//迭代求解精确的半径,确保体积分数在+-0.5%
	setROuter(ini_r, -1);
	genFusionTrussVoxel();
	double vof_current = (double)voxel.sum() / (double)voxel.rows();
	//std::cout << "迭代开始时体积分数为" << vof_current << ",半径为" << ini_r << std::endl;

	double factor = 1.0;	//步进长度因子
	int fit_iterations = 0;
	while (std::abs(vof - vof_current) > 0.005)
	{
		if (factor < 1e-10 || ++fit_iterations > 1000) return false;
		if (vof_current < vof)
		{
			//一点点增加
			while (vof - vof_current > 0.005)
			{
				if (++fit_iterations > 1000) return false;
				ini_r += 0.1 * factor;
				setROuter(ini_r, -1);
				genFusionTrussVoxel();
				vof_current = (double)voxel.sum() / (double)voxel.rows();
			}
		}
		else if (vof_current > vof)
		{
			//一点点减少
			while (vof - vof_current < -0.005)
			{
				if (++fit_iterations > 1000) return false;
				ini_r -= 0.1 * factor;
				setROuter(ini_r, -1);
				genFusionTrussVoxel();
				vof_current = (double)voxel.sum() / (double)voxel.rows();
			}
		}
		else {}
		factor = factor * 0.1;
		//std::cout << "迭代时体积分数为" << vof_current << ",半径为" << ini_r << std::endl;
	}
	std::cout << "最终体积分数为" << vof_current << ",半径为" << ini_r << std::endl;
	if (ini_r < 0.8)
		return false;
	return true;
}

void Truss::genFusionTrussVoxel()
{
	double cell_length_half = cell_length / 2;	//模型大小的一半,与一些数值参数有关

	//重建范围
	Eigen::RowVector3d VminCut = { -1 * cell_length_half, -1 * cell_length_half, -1 * cell_length_half };
	Eigen::RowVector3d VmaxCut = { cell_length_half, cell_length_half, cell_length_half };

	Eigen::RowVector3i res = (resolution * ((VmaxCut - VminCut) / (VmaxCut - VminCut).maxCoeff())).cast<int>();

	grid.resize(res(0) * res(1) * res(2), 3);

	//构建采样点的坐标
	const auto lerpCut = [&](const int di, const int d)->double
	{return VminCut(d) + (double)di / (double)(res(d) - 1) * (VmaxCut(d) - VminCut(d)); };

#pragma omp parallel for
	for (int zi = 0; zi < res(2); zi++)
	{
		const double zc = lerpCut(zi, 2);
		for (int yi = 0; yi < res(1); yi++)
		{
			const double yc = lerpCut(yi, 1);
			for (int xi = 0; xi < res(0); xi++)
			{
				const double xc = lerpCut(xi, 0);
				grid.row(xi + res(0) * (yi + res(1) * zi)) = Eigen::RowVector3d(xc, yc, zc);
			}
		}
	}

	//std::vector<Edge> edge_list = base_edges;
	std::vector<Edge> edge_list = edges;
	std::vector<Joint> joint_list = nodes;
	double smf = smooth_factor;

	//计算距离场的公式
	//看PPL的计算
	//将采样点xyz映射到子单元1/8或者1/48中,然后计算距离,只需要关照子单元中的几条杆
	const auto& function = [edge_list, joint_list, cell_length_half, smf](double x, double y, double z)->double {
		//超出边界,赋一个比较大的正值
		if (abs(x) > cell_length_half || abs(y) > cell_length_half || abs(z) > cell_length_half)
		{
			return 100;
		}
		//将采样点映射到子单元
		Eigen::Vector3d p = { x,y,z };
		getPanettaSymmetryPoint(p);

		double result_KS = 0;
		for (int i = 0; i < edge_list.size(); i++)
		{
			Eigen::Vector3d v1 = edge_list[i].v.row(0), v2 = edge_list[i].v.row(1);
			double dis = pointToSegDist(v1, v2, p) - edge_list[i].r_outer;
			result_KS += exp(-1 / smf * dis);	//smf是公式中的1/s_j
		}
		result_KS = -smf * log(result_KS);
		return result_KS;
	};
	voxel.resize(grid.rows());
	voxel.setZero();
#pragma omp parallel for
	for (int i = 0; i < grid.rows(); i++) {
		double indicator = function(grid(i, 0), grid(i, 1), grid(i, 2));
		if (indicator <= 0)
		{
			voxel(i) = 1;
		}
		else
			voxel(i) = 0;
	}
}

void Truss::genFusionTrussMesh()
{
	double cell_length_half = cell_length / 2;	//模型大小的一半,与一些数值参数有关

	//重建范围
	Eigen::RowVector3d Vmin = { -1 * cell_length_half - 0.1, -1 * cell_length_half - 0.1, -1 * cell_length_half - 0.1 };
	Eigen::RowVector3d Vmax = { cell_length_half + 0.1, cell_length_half + 0.1, cell_length_half + 0.1 };

	Eigen::RowVector3i res_mesh = (resolution_mesh * ((Vmax - Vmin) / (Vmax - Vmin).maxCoeff())).cast<int>();

	Eigen::MatrixXd GV(res_mesh(0) * res_mesh(1) * res_mesh(2), 3);

	//构建采样点的坐标
	const auto lerp = [&](const int di, const int d)->double
	{return Vmin(d) + (double)di / (double)(res_mesh(d) - 1) * (Vmax(d) - Vmin(d)); };

#pragma omp parallel for
	for (int zi = 0; zi < res_mesh(2); zi++)
	{
		const double z = lerp(zi, 2);
		for (int yi = 0; yi < res_mesh(1); yi++)
		{
			const double y = lerp(yi, 1);
			for (int xi = 0; xi < res_mesh(0); xi++)
			{
				const double x = lerp(xi, 0);
				GV.row(xi + res_mesh(0) * (yi + res_mesh(1) * zi)) = Eigen::RowVector3d(x, y, z);
			}
		}
	}

	//std::vector<Edge> edge_list = base_edges;
	std::vector<Edge> edge_list = edges;
	std::vector<Joint> joint_list = nodes;
	double smf = smooth_factor;

	//计算距离场的公式
	//看PPL的计算
	//将采样点xyz映射到子单元1/8或者1/48中,然后计算距离,只需要关照子单元中的几条杆
	const auto& function = [edge_list, joint_list, cell_length_half, smf](double x, double y, double z)->double {
		//超出边界,赋一个比较大的正值
		if (abs(x) > cell_length_half || abs(y) > cell_length_half || abs(z) > cell_length_half)
		{
			return 100;
		}
		//将采样点映射到子单元
		Eigen::Vector3d p = { x,y,z };
		getPanettaSymmetryPoint(p);

		double result_KS = 0;
		for (int i = 0; i < edge_list.size(); i++)
		{
			Eigen::Vector3d v1 = edge_list[i].v.row(0), v2 = edge_list[i].v.row(1);
			double dis = pointToSegDist(v1, v2, p) - edge_list[i].r_outer;
			result_KS += exp(-1 / smf * dis);	//smf是公式中的1/s_j
		}
		result_KS = -smf * log(result_KS);
		return result_KS;
	};

	Eigen::VectorXd B(GV.rows());
#pragma omp parallel for
	for (int i = 0; i < GV.rows(); i++) {
		double indicator = function(GV(i, 0), GV(i, 1), GV(i, 2));
		B(i) = indicator;
	}
	igl::copyleft::marching_cubes(B, GV, res_mesh(0), res_mesh(1), res_mesh(2), V, F);

	//平移中心到原点
	if (V.rows() == 0 || !V.allFinite())
		throw std::runtime_error("No finite isosurface; check the full-cell skeleton, radius and mesh resolution");
	Eigen::Vector3d minV = V.colwise().minCoeff();
	Eigen::Vector3d maxV = V.colwise().maxCoeff();

	Eigen::Vector3d move = minV + (maxV - minV) / 2.0;
	for (int j = 0; j < V.rows(); j++)
	{
		V(j, 0) -= move(0);
		V(j, 1) -= move(1);
		V(j, 2) -= move(2);
	}

	//缩放为[-1,-1,-1]到[1,1,1]的边长为2的立方体
	Eigen::Vector3d factor = (maxV - minV) / 2.0;
	double xfac = (1.0 / factor(0));
	double yfac = (1.0 / factor(1));
	double zfac = (1.0 / factor(2));
	for (int j = 0; j < V.rows(); j++)
	{

		V(j, 0) *= xfac;
		V(j, 1) *= yfac;
		V(j, 2) *= zfac;
	}
}


void Truss::genTilingTrussMesh(int size_x, int size_y, int size_z)
{
	//应注意重建范围 [0,0,0]->[cell_length*size_x,cell_length*size_y,cell_length*size_z] 再加入边界一小圈
	Eigen::RowVector3d Vmin = { 0.0, 0.0, 0.0 };
	Eigen::RowVector3d Vmax = { cell_length * size_x, cell_length * size_y, cell_length * size_z };

	//这里减掉一些是因为每两个格子会重复用一个采样点
	Eigen::RowVector3i res_mesh = { resolution_mesh * size_x - (size_x - 1),resolution_mesh * size_y - (size_y - 1),resolution_mesh * size_z - (size_z - 1) };


	
	Eigen::MatrixXd GV((res_mesh(0) + 2) * (res_mesh(1) + 2) * (res_mesh(2) + 2), 3);	//算上边界

	//构建采样点的坐标
	const auto lerp = [&](const int di, const int d)->double
	{return Vmin(d) + (double)di / (double)(res_mesh(d) - 1) * (Vmax(d) - Vmin(d)); };

	double boundary_length = 0.1;

#pragma omp parallel for
	for (int zi = 0; zi < res_mesh(2) + 2; zi++)
	{
		double z = 0;
		if (zi == 0)
			z = Vmin(2) - boundary_length;
		else if (zi == res_mesh(2) + 1)
			z = Vmax(2) + boundary_length;
		else
			z = lerp(zi, 2);
		for (int yi = 0; yi < res_mesh(1) + 2; yi++)
		{
			double y = 0;
			if (yi == 0)
				y = Vmin(1) - boundary_length;
			else if (yi == res_mesh(1) + 1)
				y = Vmax(1) + boundary_length;
			else
				y = lerp(yi, 1);
			for (int xi = 0; xi < res_mesh(0) + 2; xi++)
			{
				double x = 0;
				if (xi == 0)
					x = Vmin(0) - boundary_length;
				else if (xi == res_mesh(0) + 1)
					x = Vmax(0) + boundary_length;
				else
					x = lerp(xi, 0);
				GV.row(xi + (res_mesh(0) + 2) * (yi + (res_mesh(1) + 2) * zi)) = Eigen::RowVector3d(x, y, z);
			}
		}
	}
	//接下来构建边界
	/*
	int GV_index = res_mesh(0) * res_mesh(1) * res_mesh(2);
	//std::cout << "GV_size=" << GV.rows() << ",GV_index=" << GV_index << std::endl;
	//顶部和底部 2xy
	for (int xi = 0; xi < res_mesh(0); xi++)
	{
		const double x = lerp(xi, 0);
		for (int yi = 0; yi < res_mesh(1); yi++)
		{
			const double y = lerp(yi, 1);
			GV.row(GV_index) = Eigen::RowVector3d(x, y, Vmin(2) - boundary_length);
			GV_index++;
			GV.row(GV_index) = Eigen::RowVector3d(x, y, Vmax(2) + boundary_length);
			GV_index++;
		}
	}
	//正面和背面 2xz
	for (int xi = 0; xi < res_mesh(0); xi++)
	{
		const double x = lerp(xi, 0);
		for (int zi = 0; zi < res_mesh(2); zi++)
		{
			const double z = lerp(zi, 2);
			GV.row(GV_index) = Eigen::RowVector3d(x, Vmin(1) - boundary_length, z);
			GV_index++;
			GV.row(GV_index) = Eigen::RowVector3d(x, Vmax(1) + boundary_length, z);
			GV_index++;
		}
	}
	//左右两侧 2yz
	for (int yi = 0; yi < res_mesh(1); yi++)
	{
		const double y = lerp(yi, 1);
		for (int zi = 0; zi < res_mesh(2); zi++)
		{
			const double z = lerp(zi, 2);
			GV.row(GV_index) = Eigen::RowVector3d(Vmin(0) - boundary_length, y, z);
			GV_index++;
			GV.row(GV_index) = Eigen::RowVector3d(Vmax(0) + boundary_length, y, z);
			GV_index++;
		}
	}
	//std::cout << "GV_size=" << GV.rows() << ",GV_index=" << GV_index << std::endl;
	//12条边
	for (int xi = 0; xi < res_mesh(0); xi++)
	{
		const double x = lerp(xi, 0);
		GV.row(GV_index) = Eigen::RowVector3d(x, Vmin(1) - boundary_length, Vmin(2) - boundary_length);
		GV_index++;
		GV.row(GV_index) = Eigen::RowVector3d(x, Vmin(1) - boundary_length, Vmax(2) + boundary_length);
		GV_index++;
		GV.row(GV_index) = Eigen::RowVector3d(x, Vmax(1) + boundary_length, Vmin(2) - boundary_length);
		GV_index++;
		GV.row(GV_index) = Eigen::RowVector3d(x, Vmax(1) + boundary_length, Vmax(2) + boundary_length);
		GV_index++;
	}
	for (int yi = 0; yi < res_mesh(1); yi++)
	{
		const double y = lerp(yi, 1);
		GV.row(GV_index) = Eigen::RowVector3d(Vmin(0) - boundary_length, y, Vmin(2) - boundary_length);
		GV_index++;
		GV.row(GV_index) = Eigen::RowVector3d(Vmin(0) - boundary_length, y, Vmax(2) + boundary_length);
		GV_index++;
		GV.row(GV_index) = Eigen::RowVector3d(Vmax(0) + boundary_length, y, Vmin(2) - boundary_length);
		GV_index++;
		GV.row(GV_index) = Eigen::RowVector3d(Vmax(0) + boundary_length, y, Vmax(2) + boundary_length);
		GV_index++;
	}
	for (int zi = 0; zi < res_mesh(2); zi++)
	{
		const double z = lerp(zi, 0);
		GV.row(GV_index) = Eigen::RowVector3d(Vmin(0) - boundary_length, Vmin(1) - boundary_length, z);
		GV_index++;
		GV.row(GV_index) = Eigen::RowVector3d(Vmin(0) - boundary_length, Vmax(1) + boundary_length, z);
		GV_index++;
		GV.row(GV_index) = Eigen::RowVector3d(Vmax(0) + boundary_length, Vmin(1) - boundary_length, z);
		GV_index++;
		GV.row(GV_index) = Eigen::RowVector3d(Vmax(0) + boundary_length, Vmax(1) + boundary_length, z);
		GV_index++;
	}
	//std::cout << "GV_size=" << GV.rows() << ",GV_index=" << GV_index << std::endl;
	//8个顶点
	GV.row(GV_index) = Eigen::RowVector3d(Vmin(0) - boundary_length, Vmin(1) - boundary_length, Vmin(2) - boundary_length);
	GV_index++;
	GV.row(GV_index) = Eigen::RowVector3d(Vmin(0) - boundary_length, Vmin(1) - boundary_length, Vmax(2) + boundary_length);
	GV_index++;
	GV.row(GV_index) = Eigen::RowVector3d(Vmin(0) - boundary_length, Vmax(1) + boundary_length, Vmin(2) - boundary_length);
	GV_index++;
	GV.row(GV_index) = Eigen::RowVector3d(Vmin(0) - boundary_length, Vmax(1) + boundary_length, Vmax(2) + boundary_length);
	GV_index++;
	GV.row(GV_index) = Eigen::RowVector3d(Vmax(0) + boundary_length, Vmin(1) - boundary_length, Vmin(2) - boundary_length);
	GV_index++;
	GV.row(GV_index) = Eigen::RowVector3d(Vmax(0) + boundary_length, Vmin(1) - boundary_length, Vmax(2) + boundary_length);
	GV_index++;
	GV.row(GV_index) = Eigen::RowVector3d(Vmax(0) + boundary_length, Vmax(1) + boundary_length, Vmin(2) - boundary_length);
	GV_index++;
	GV.row(GV_index) = Eigen::RowVector3d(Vmax(0) + boundary_length, Vmax(1) + boundary_length, Vmax(2) + boundary_length);
	GV_index++;

	std::cout << "GV_size=" << GV.rows() << ",GV_index=" << GV_index << std::endl;*/

	//std::vector<Edge> edge_list = base_edges;
	std::vector<Edge> edge_list = edges;
	std::vector<Joint> joint_list = nodes;
	double cell_len = cell_length;
	//double smf = smooth_factor;

	//计算距离场的公式
	//看PPL的计算
	//将采样点xyz映射到子单元1/8或者1/48中,然后计算距离,只需要关照子单元中的几条杆
	const auto& function = [edge_list, joint_list, Vmin, Vmax, cell_len](double x, double y, double z)->double {
		//超出边界,赋一个比较大的正值
		if (x<Vmin(0) || x>Vmax(0) || y<Vmin(1) || y>Vmax(1) || z<Vmin(2) || z>Vmax(2))
		{
			return 100;
		}
		//将采样点映射到子单元,这里注意先取余到第一个单元[0,0,0]->[cell_length,cell_length,cell_length],再平移到以0为原点
		while (x > cell_len)
			x = x - cell_len;
		while (y > cell_len)
			y = y - cell_len;
		while (z > cell_len)
			z = z - cell_len;
		Eigen::Vector3d p = { x - cell_len / 2.0,y - cell_len / 2.0,z - cell_len / 2.0 };
		getPanettaSymmetryPoint(p);

		std::vector<double> result;	//存放每个joint的KS函数值,到最后再套一层总的KS函数

		for (int i = 0; i < joint_list.size(); i++)
		{
			double result_joint = 0;		//当前joint的KS函数值
			double smf = joint_list[i].smooth_factor;
			for (int j = 0; j < joint_list[i].associated_edges.size(); j++)
			{
				int index = joint_list[i].associated_edges[j];	//相邻边在edge_list中的索引
				Eigen::Vector3d v1 = edge_list[index].v.row(0), v2 = edge_list[index].v.row(1);
				double dis = pointToSegDist(v1, v2, p) - edge_list[index].r_outer;
				result_joint += exp(-1 / smf * dis);	//smf是公式中的1/s_j
			}
			result_joint = -smf * log(result_joint);
			result.push_back(result_joint);
		}

		double sum_KS = 0;	//总的KS函数值
		double sum_smf = 0.5;
		for (int i = 0; i < result.size(); i++)
		{
			sum_KS += exp(-1 / sum_smf * result[i]);
		}
		sum_KS = -sum_smf * log(sum_KS);
		return sum_KS;
	};

	Eigen::VectorXd B(GV.rows());
#pragma omp parallel for
	for (int i = 0; i < GV.rows(); i++) {
		double indicator = function(GV(i, 0), GV(i, 1), GV(i, 2));
		B(i) = indicator;
	}
	igl::copyleft::marching_cubes(B, GV, (res_mesh(0) + 2), (res_mesh(1) + 2), (res_mesh(2) + 2), V, F);

	//平移min到原点
	if (V.rows() == 0 || !V.allFinite())
		throw std::runtime_error("No finite isosurface; check the full-cell skeleton, radius and mesh resolution");
	Eigen::Vector3d minV = V.colwise().minCoeff();
	Eigen::Vector3d maxV = V.colwise().maxCoeff();

	//Eigen::Vector3d move = minV + (maxV - minV) / 2.0;
	Eigen::Vector3d move = minV;
	for (int j = 0; j < V.rows(); j++)
	{
		V(j, 0) -= move(0);
		V(j, 1) -= move(1);
		V(j, 2) -= move(2);
	}

	//缩放为[0,0,0]开始的边长为1的立方体
	Eigen::Vector3d factor = maxV - minV;
	double xfac = ((double)size_x / factor(0));
	double yfac = ((double)size_y / factor(1));
	double zfac = ((double)size_z / factor(2));
	for (int j = 0; j < V.rows(); j++)
	{

		V(j, 0) *= xfac;
		V(j, 1) *= yfac;
		V(j, 2) *= zfac;
	}
}





//判断E是否在同一个平面,如果是,点p是否在面上
//E:都是Edges的向量,末尾点-起始点 3*n矩阵
//p:3*1向量,p也要减去起始点
void checkCoplanar(Eigen::MatrixXd E, Eigen::Vector3d p, bool& isCoplanar, bool& isInPlane)
{
	int n = E.cols();
	if (n <= 1)
	{
		isCoplanar = true;
		isInPlane = false;
		return;
	}
	isCoplanar = true;
	isInPlane = false;

	//选择前两个向量来确定平面
	Eigen::Vector3d e1 = E.col(0);
	Eigen::Vector3d e2 = E.col(1);
	Eigen::Vector3d normal = e1.cross(e2);

	//检查剩余向量是否共面
	for (int i = 2; i < n; ++i) {
		Eigen::Vector3d ei = E.col(i);
		if (abs(ei.dot(normal)) > 1e-4) {  // 允许浮点误差
			isCoplanar = false;
			return;
		}
	}
	//如果所有向量共面,检查点p是否在平面内
	if (isCoplanar) {
		double d = abs(p.dot(normal)) / normal.norm();
		if (d < 1e-5) {  // 距离接近零，考虑浮点数误差
			isInPlane = true;
		}
	}
}


//判断点p是否在E的内侧,即存在p=a1e1+a2e2+...+anen,ai>=0
//E:都是Edges的向量,末尾点-起始点 3*n矩阵
//p:3*1向量,p也要减去起始点
//出现多解时,可能出错,所以用nnls计算非负解,然后看这个解可不可行
bool checkNonNegativeCombination(Eigen::MatrixXd E, Eigen::Vector3d p)
{
	//先对E进行标准化
	for (int i = 0; i < E.cols(); i++)
	{
		E.col(i) = E.col(i).normalized();
	}
	p = p.normalized();
	//p = p * 100;

	//std::cout << E << std::endl;
	//std::cout << p << std::endl;


	//特殊处理:所有边在一个平面上
	bool isCoplanar = false, isInPlane = false;
	checkCoplanar(E, p, isCoplanar, isInPlane);
	if (isCoplanar && !isInPlane)	//不在平面内,肯定无法组合
		return false;

	//检查是否存在完全相反的e1和e2,有则不用判断a1a2符号
	int n = E.cols();
	std::vector<int> needCheck;
	for (int i = 0; i < n; i++)
		needCheck.push_back(1);

	for (int i = 0; i < n - 1; ++i) {
		for (int j = i + 1; j < n; ++j) {
			Eigen::Vector3d ei = E.col(i);
			Eigen::Vector3d ej = E.col(j);

			// 需要叉乘为0且点乘小于0
			if (ei.cross(ej).norm() == 0 && ei.dot(ej) < 0) {
				needCheck[i] = 0;
				needCheck[j] = 0;
			}
		}
	}
	

	Eigen::NNLS<Eigen::MatrixXd> nnls(E);
	nnls.solve(p);
	Eigen::VectorXd a = nnls.x();

	Eigen::VectorXd r = p - E * a;
	double residualNorm = r.norm();
	if (residualNorm > 1e-3)
		return false;

	//求解线性方程组E*a=p
	//Eigen::VectorXd a = E.colPivHouseholderQr().solve(p);

	//如果有负的,先试一下能不能转为其他基向量

	//检查解向量是否全部大于等于0
	for (int i = 0; i < n; i++)
	{
		if (needCheck[i] == 1 && a(i) < 0)
		{
			return false;
		}
	}
	return true;
}

//点p到射线oa的距离
double distancePointToRay(Eigen::Vector3d o, Eigen::Vector3d a, Eigen::Vector3d p)
{
	// 计算向量 oa 和 op
	Eigen::Vector3d oa = a - o;
	Eigen::Vector3d op = p - o;

	// 计算向量 oa 的长度
	double oa_length = oa.norm();

	// 如果 oa 的长度为零，射线退化成点
	if (oa_length == 0) {
		std::cout << "点 o 和点 a 重合，射线退化为点。" << std::endl;
	}

	// 计算单位向量 u
	Eigen::Vector3d u = oa / oa_length;

	// 计算向量 op 在向量 oa 方向上的投影
	double op_dot_u = op.dot(u);

	double distance;
	if (op_dot_u <= 0) {
		// 如果投影小于等于 0，点 p 在点 o 的左侧
		distance = op.norm();
	}
	else {
		// 计算投影向量
		Eigen::Vector3d op_parallel_oa = op_dot_u * u;

		// 计算垂直于 oa 的向量
		Eigen::Vector3d op_perp_oa = op - op_parallel_oa;

		// 计算距离
		distance = op_perp_oa.norm();
	}

	return distance;
}

bool pointInTriangle(const Eigen::Vector3d& p, const Eigen::Vector3d& o, const Eigen::Vector3d& a, const Eigen::Vector3d& b) {
	Eigen::Vector3d v0 = b - o;
	Eigen::Vector3d v1 = a - o;
	Eigen::Vector3d v2 = p - o;

	double dot00 = v0.dot(v0);
	double dot01 = v0.dot(v1);
	double dot02 = v0.dot(v2);
	double dot11 = v1.dot(v1);
	double dot12 = v1.dot(v2);

	double invDenom = 1.0 / (dot00 * dot11 - dot01 * dot01);
	double u = (dot11 * dot02 - dot01 * dot12) * invDenom;
	double v = (dot00 * dot12 - dot01 * dot02) * invDenom;

	return (u >= 0) && (v >= 0) && (u + v < 1);
}


//distance1:点p到两条射线oa和ob组成的部分平面的距离
//distance2:点p到三角形oab的距离
std::pair<double, double> pointToPartialPlaneDistance(Eigen::Vector3d p, Eigen::Vector3d o, Eigen::Vector3d a, Eigen::Vector3d b)
{
	// 向量
	Eigen::Vector3d oa = a - o;
	Eigen::Vector3d ob = b - o;
	Eigen::Vector3d op = p - o;

	// 计算法向量
	Eigen::Vector3d normal = oa.cross(ob);
	double normal_length = normal.norm();

	// 投影点到平面的距离
	double distance_to_plane = std::abs(normal.dot(op)) / normal_length;

	// 投影点
	Eigen::Vector3d p_proj = p - (normal.dot(op) / normal_length / normal_length) * normal;

	// 向量
	Eigen::Vector3d op_proj = p_proj - o;
	op_proj.normalize();

	// 检查投影点是否在扇形平面内
	// 首先p和b在oa同一侧
	double dis1 = 0;
	bool dis1_flag = false;
	if (oa.cross(op_proj).dot(normal) >= 0) {
		double angle_AOB = std::acos(oa.normalized().dot(ob.normalized()));
		double angle_AOP = std::acos(oa.normalized().dot(op_proj));
		if (angle_AOP <= angle_AOB) {
			dis1 = distance_to_plane;
			dis1_flag = true;
		}
	}

	// 计算到射线 oa 和 ob 的距离
	if (!dis1_flag)
	{
		double distance_to_oa = distancePointToRay(o, a, p);
		double distance_to_ob = distancePointToRay(o, b, p);
		dis1 = std::min(distance_to_oa, distance_to_ob);
	}
	
	//检查投影点是否在三角形内
	double dis2 = 0;
	bool dis2_flag = false;
	if (pointInTriangle(p_proj, o, a, b))
	{
		dis2 = distance_to_plane;
		dis2_flag = true;
	}
	if (!dis2_flag)
	{
		double distance_to_oa = pointToSegDist(o, a, p);
		double distance_to_ob = pointToSegDist(o, b, p);
		double distance_to_ab = pointToSegDist(a, b, p);
		dis2 = std::min(distance_to_oa, distance_to_ob);
		dis2 = std::min(dis2, distance_to_ab);
	}
	return std::pair<double, double>(dis1, dis2);
}

//distance:点p到三角形oab的距离
double pointToTriangleDistance(Eigen::Vector3d p, Eigen::Vector3d o, Eigen::Vector3d a, Eigen::Vector3d b)
{
	// 向量
	Eigen::Vector3d oa = a - o;
	Eigen::Vector3d ob = b - o;
	Eigen::Vector3d op = p - o;

	// 计算法向量
	Eigen::Vector3d normal = oa.cross(ob);
	double normal_length = normal.norm();

	// 投影点到平面的距离
	double distance_to_plane = std::abs(normal.dot(op)) / normal_length;

	// 投影点
	Eigen::Vector3d p_proj = p - (normal.dot(op) / normal_length / normal_length) * normal;

	//检查投影点是否在三角形内
	double dis = 0;
	if (pointInTriangle(p_proj, o, a, b))
	{
		dis = distance_to_plane;
		return dis;
	}
	//否则计算到三条线段的距离
	double distance_to_oa = pointToSegDist(o, a, p);
	double distance_to_ob = pointToSegDist(o, b, p);
	double distance_to_ab = pointToSegDist(a, b, p);
	dis = std::min(distance_to_oa, distance_to_ob);
	dis = std::min(dis, distance_to_ab);
	return dis;
}




bool Truss::setSolidVolumeFractionModified(double vof)
{
	//计算当前体积分数
	genSolidTrussVoxelModified();
	double vof_1 = (double)voxel.sum() / (double)voxel.rows();
	std::cout << "初始体积分数为" << vof_1 << ",半径为" << r_outer << std::endl;

	//将d_outer加一点点,然后计算斜率
	setROuter(r_outer + 0.5, -1);
	genSolidTrussVoxelModified();
	double vof_2 = (double)voxel.sum() / (double)voxel.rows();
	std::cout << "初始体积分数为" << vof_2 << ",半径为" << r_outer << std::endl;

	double k = ((vof_2 - vof_1) * 1.0) / 0.5;//斜率
	if (!std::isfinite(k)) return false;
	double ini_r = std::abs(k) < 1e-14 ? r_outer : ((vof - vof_2) / k + r_outer);//粗略计算目标相对密度的半径
	if (ini_r < 0)
		ini_r = 1.0;		//从最小开始

	//迭代求解精确的半径,确保体积分数在+-0.5%
	setROuter(ini_r, -1);
	genSolidTrussVoxelModified();
	double vof_current = (double)voxel.sum() / (double)voxel.rows();
	std::cout << "迭代开始时体积分数为" << vof_current << ",半径为" << ini_r << std::endl;

	double factor = 1.0;	//步进长度因子
	int fit_iterations = 0;
	while (std::abs(vof - vof_current) > 0.005)
	{
		if (factor < 1e-10 || ++fit_iterations > 1000) return false;
		if (vof_current < vof)
		{
			//一点点增加
			while (vof - vof_current > 0.005)
			{
				if (++fit_iterations > 1000) return false;
				ini_r += 0.5 * factor;
				setROuter(ini_r, -1);
				genSolidTrussVoxelModified();
				vof_current = (double)voxel.sum() / (double)voxel.rows();
			}
		}
		else if (vof_current > vof)
		{
			//一点点减少
			while (vof - vof_current < -0.005)
			{
				if (++fit_iterations > 1000) return false;
				ini_r -= 0.5 * factor;
				setROuter(ini_r, -1);
				genSolidTrussVoxelModified();
				vof_current = (double)voxel.sum() / (double)voxel.rows();
			}
		}
		else {}
		factor = factor * 0.1;
		std::cout << "迭代时体积分数为" << vof_current << ",半径为" << ini_r << std::endl;
	}
	std::cout << "最终体积分数为" << vof_current << ",半径为" << ini_r << std::endl;
	if (ini_r < 0.8)
		return false;
	return true;
}


void Truss::genSolidTrussVoxelModified()
{
	double cell_length_half = cell_length / 2;	//模型大小的一半,与一些数值参数有关

	//重建范围
	Eigen::RowVector3d VminCut = { -1 * cell_length_half, -1 * cell_length_half, -1 * cell_length_half };
	Eigen::RowVector3d VmaxCut = { cell_length_half, cell_length_half, cell_length_half };

	Eigen::RowVector3i res = (resolution * ((VmaxCut - VminCut) / (VmaxCut - VminCut).maxCoeff())).cast<int>();

	grid.resize(res(0) * res(1) * res(2), 3);

	//构建采样点的坐标
	const auto lerpCut = [&](const int di, const int d)->double
	{return VminCut(d) + (double)di / (double)(res(d) - 1) * (VmaxCut(d) - VminCut(d)); };

#pragma omp parallel for
	for (int zi = 0; zi < res(2); zi++)
	{
		const double zc = lerpCut(zi, 2);
		for (int yi = 0; yi < res(1); yi++)
		{
			const double yc = lerpCut(yi, 1);
			for (int xi = 0; xi < res(0); xi++)
			{
				const double xc = lerpCut(xi, 0);
				grid.row(xi + res(0) * (yi + res(1) * zi)) = Eigen::RowVector3d(xc, yc, zc);
			}
		}
	}

	//std::vector<Edge> edge_list = base_edges;
	std::vector<Edge> edge_list = edges;
	std::vector<Joint> joint_list = nodes;
	//double smf = smooth_factor;

	//计算距离场的公式
	//看PPL的计算
	//将采样点xyz映射到子单元1/8或者1/48中,然后计算距离,只需要关照子单元中的几条杆
	const auto& function = [edge_list, joint_list, cell_length_half](double x, double y, double z)->double {
		//超出边界,赋一个比较大的正值
		if (abs(x) > cell_length_half || abs(y) > cell_length_half || abs(z) > cell_length_half)
		{
			return 100;
		}
		//将采样点映射到子单元
		Eigen::Vector3d p = { x,y,z };
		getPanettaSymmetryPoint(p);

		std::vector<double> result;	//存放每个joint的KS函数值,到最后再套一层总的KS函数

		//std::cout << joint_list.size() << std::endl;
		for (int i = 0; i < joint_list.size(); i++)
		{
			double result_joint = 0;		//当前joint的KS函数值
			double smf = joint_list[i].smooth_factor;
			for (int j = 0; j < joint_list[i].associated_edges.size(); j++)
			{
				int index = joint_list[i].associated_edges[j];	//相邻边在edge_list中的索引
				Eigen::Vector3d v1 = edge_list[index].v.row(0), v2 = edge_list[index].v.row(1);
				double dis = pointToSegDist(v1, v2, p) - edge_list[index].r_outer;
				result_joint += exp(-1 / smf * dis);	//smf是公式中的1/s_j
			}
			result_joint = -smf * log(result_joint);

			//判断该joint的n条边是否在空间同一侧,如果是则继续,否则正常计算KS函数值
			if (joint_list[i].is_same_side)
			{
				//判断点p是否在内侧
				Eigen::Vector3d pp = p - joint_list[i].center;
				int n = joint_list[i].associated_edges.size();
				Eigen::MatrixXd E(3, n);

				for (int j = 0; j < n; j++)
				{
					int index = joint_list[i].associated_edges[j];	//相邻边在edge_list中的索引
					Eigen::Vector3d v1 = edge_list[index].v.row(0), v2 = edge_list[index].v.row(1);
					if (isEquation(joint_list[i].center, v1))
						E.col(j) = v2 - joint_list[i].center;
					else if (isEquation(joint_list[i].center, v2))
						E.col(j) = v1 - joint_list[i].center;
					else
						std::cout << "Error in getTopology() calculate is_same_side" << std::endl;
				}
				bool isIn=checkNonNegativeCombination(E, pp);
				//如果在外侧,需要切割
				if (!isIn)
				{
					//计算点p对应哪个边界,即计算到所有扇形面距离减去r,取min
					double minDis = 1e9;
					for (int j = 0; j < n - 1; j++)
					{
						Eigen::Vector3d ej;
						int index_j = joint_list[i].associated_edges[j];	//相邻边在edge_list中的索引
						Eigen::Vector3d v1 = edge_list[index_j].v.row(0), v2 = edge_list[index_j].v.row(1);
						if (isEquation(joint_list[i].center, v1))
							ej = v2;
						else if (isEquation(joint_list[i].center, v2))
							ej = v1;
						else
							std::cout << "Error in getTopology() calculate is_same_side" << std::endl;
						for (int k = j + 1; k < n; k++)
						{
							Eigen::Vector3d ek;
							int index_k = joint_list[i].associated_edges[k];	//相邻边在edge_list中的索引
							Eigen::Vector3d v1 = edge_list[index_k].v.row(0), v2 = edge_list[index_k].v.row(1);
							if (isEquation(joint_list[i].center, v1))
								ek = v2;
							else if (isEquation(joint_list[i].center, v2))
								ek = v1;
							else
								std::cout << "Error in getTopology() calculate is_same_side" << std::endl;
							
							/*double dis = pointToTriangleDistance(p, joint_list[i].center, ej, ek);
							dis = dis - std::max(edge_list[index_j].r_outer, edge_list[index_k].r_outer);
							if (dis < minDis)
								minDis = dis;*/

							std::pair<double, double> dis_pair = pointToPartialPlaneDistance(p, joint_list[i].center, ej, ek);
							double dis1 = dis_pair.first - std::max(edge_list[index_j].r_outer, edge_list[index_k].r_outer);
							double dis2 = dis_pair.second - std::max(edge_list[index_j].r_outer, edge_list[index_k].r_outer);
							if (dis1 < minDis)
								minDis = dis1;
							//if (dis2 < minDis2)
							//	minDis2 = dis2;
						}
					}

					if (minDis >= 0)
					{
						result_joint = std::max(result_joint, minDis);
					}
					else if (result_joint < 0)
						result_joint = result_joint * (-1) * tanh(minDis);

				}
			}
			result.push_back(result_joint);
		}

		double sum_KS = 0;	//总的KS函数值
		double sum_smf = 0.1;
		for (int i = 0; i < result.size(); i++)
		{
			sum_KS += exp(-1 / sum_smf * result[i]);
		}
		sum_KS = -sum_smf * log(sum_KS);
		return sum_KS;
	};
	voxel.resize(grid.rows());
	voxel.setZero();
#pragma omp parallel for
	for (int i = 0; i < grid.rows(); i++) {
		double indicator = function(grid(i, 0), grid(i, 1), grid(i, 2));
		if (indicator <= 0)
		{
			voxel(i) = 1;
		}
		else
			voxel(i) = 0;
	}
}


void Truss::genSolidTrussMeshModified()
{
	double cell_length_half = cell_length / 2;	//模型大小的一半,与一些数值参数有关

	//重建范围
	Eigen::RowVector3d Vmin = { -1 * cell_length_half - 0.1, -1 * cell_length_half - 0.1, -1 * cell_length_half - 0.1 };
	Eigen::RowVector3d Vmax = { cell_length_half + 0.1, cell_length_half + 0.1, cell_length_half + 0.1 };

	Eigen::RowVector3i res_mesh = (resolution_mesh * ((Vmax - Vmin) / (Vmax - Vmin).maxCoeff())).cast<int>();

	Eigen::MatrixXd GV(res_mesh(0) * res_mesh(1) * res_mesh(2), 3);

	//构建采样点的坐标
	const auto lerp = [&](const int di, const int d)->double
	{return Vmin(d) + (double)di / (double)(res_mesh(d) - 1) * (Vmax(d) - Vmin(d)); };

#pragma omp parallel for
	for (int zi = 0; zi < res_mesh(2); zi++)
	{
		const double z = lerp(zi, 2);
		for (int yi = 0; yi < res_mesh(1); yi++)
		{
			const double y = lerp(yi, 1);
			for (int xi = 0; xi < res_mesh(0); xi++)
			{
				const double x = lerp(xi, 0);
				GV.row(xi + res_mesh(0) * (yi + res_mesh(1) * zi)) = Eigen::RowVector3d(x, y, z);
			}
		}
	}

	//std::vector<Edge> edge_list = base_edges;
	std::vector<Edge> edge_list = edges;
	std::vector<Joint> joint_list = nodes;
	//double smf = smooth_factor;

	//计算距离场的公式
	//看PPL的计算
	//将采样点xyz映射到子单元1/8或者1/48中,然后计算距离,只需要关照子单元中的几条杆
	const auto& function = [edge_list, joint_list, cell_length_half](double x, double y, double z)->double {
		//超出边界,赋一个比较大的正值
		if (abs(x) > cell_length_half || abs(y) > cell_length_half || abs(z) > cell_length_half)
		{
			return 100;
		}
		//将采样点映射到子单元
		Eigen::Vector3d p = { x,y,z };
		getPanettaSymmetryPoint(p);

		std::vector<double> result;	//存放每个joint的KS函数值,到最后再套一层总的KS函数

		for (int i = 0; i < joint_list.size(); i++)
		{
			double result_joint = 0;		//当前joint的KS函数值
			double smf = joint_list[i].smooth_factor;
			for (int j = 0; j < joint_list[i].associated_edges.size(); j++)
			{
				int index = joint_list[i].associated_edges[j];	//相邻边在edge_list中的索引
				Eigen::Vector3d v1 = edge_list[index].v.row(0), v2 = edge_list[index].v.row(1);
				double dis = pointToSegDist(v1, v2, p) - edge_list[index].r_outer;
				result_joint += exp(-1 / smf * dis);	//smf是公式中的1/s_j
			}
			result_joint = -smf * log(result_joint);

			//判断该joint的n条边是否在空间同一侧,如果是则继续,否则正常计算KS函数值
			if (joint_list[i].is_same_side)
			{
				//判断点p是否在内侧
				Eigen::Vector3d pp = p - joint_list[i].center;
				int n = joint_list[i].associated_edges.size();
				Eigen::MatrixXd E(3, n);

				for (int j = 0; j < n; j++)
				{
					int index = joint_list[i].associated_edges[j];	//相邻边在edge_list中的索引
					Eigen::Vector3d v1 = edge_list[index].v.row(0), v2 = edge_list[index].v.row(1);
					if (isEquation(joint_list[i].center, v1))
						E.col(j) = v2 - joint_list[i].center;
					else if (isEquation(joint_list[i].center, v2))
						E.col(j) = v1 - joint_list[i].center;
					else
						std::cout << "Error in getTopology() calculate is_same_side" << std::endl;
				}
				bool isIn = checkNonNegativeCombination(E, pp);
				//如果在外侧,需要切割
				if (!isIn)
				{
					//计算点p对应哪个边界,即计算到所有扇形面距离减去r,取min
					double minDis = 1e9;
					for (int j = 0; j < n - 1; j++)
					{
						Eigen::Vector3d ej;
						int index_j = joint_list[i].associated_edges[j];	//相邻边在edge_list中的索引
						Eigen::Vector3d v1 = edge_list[index_j].v.row(0), v2 = edge_list[index_j].v.row(1);
						if (isEquation(joint_list[i].center, v1))
							ej = v2;
						else if (isEquation(joint_list[i].center, v2))
							ej = v1;
						else
							std::cout << "Error in getTopology() calculate is_same_side" << std::endl;
						for (int k = j + 1; k < n; k++)
						{
							Eigen::Vector3d ek;
							int index_k = joint_list[i].associated_edges[k];	//相邻边在edge_list中的索引
							Eigen::Vector3d v1 = edge_list[index_k].v.row(0), v2 = edge_list[index_k].v.row(1);
							if (isEquation(joint_list[i].center, v1))
								ek = v2;
							else if (isEquation(joint_list[i].center, v2))
								ek = v1;
							else
								std::cout << "Error in getTopology() calculate is_same_side" << std::endl;

							/*double dis = pointToTriangleDistance(p, joint_list[i].center, ej, ek);
							dis = dis - std::max(edge_list[index_j].r_outer, edge_list[index_k].r_outer);
							if (dis < minDis)
								minDis = dis;*/

							std::pair<double, double> dis_pair = pointToPartialPlaneDistance(p, joint_list[i].center, ej, ek);
							double dis1 = dis_pair.first - std::max(edge_list[index_j].r_outer, edge_list[index_k].r_outer);
							double dis2 = dis_pair.second - std::max(edge_list[index_j].r_outer, edge_list[index_k].r_outer);
							if (dis1 < minDis)
								minDis = dis1;
							//if (dis2 < minDis2)
							//	minDis2 = dis2;
						}
					}

					if (minDis >= 0)
					{
						result_joint = std::max(result_joint, minDis);
					}
					else if (result_joint < 0)
						result_joint = result_joint * (-1) * tanh(minDis);

				}
			}
			result.push_back(result_joint);
		}

		double sum_KS = 0;	//总的KS函数值
		double sum_smf = 0.1;
		for (int i = 0; i < result.size(); i++)
		{
			sum_KS += exp(-1 / sum_smf * result[i]);
		}
		sum_KS = -sum_smf * log(sum_KS);
		return sum_KS;
	};

	Eigen::VectorXd B(GV.rows());
#pragma omp parallel for
	for (int i = 0; i < GV.rows(); i++) {
		double indicator = function(GV(i, 0), GV(i, 1), GV(i, 2));
		B(i) = indicator;
	}

	//输出一下距离场B
	// SDF values can be exported here for debugging.



	igl::copyleft::marching_cubes(B, GV, res_mesh(0), res_mesh(1), res_mesh(2), V, F);

	//平移中心到原点
	if (V.rows() == 0 || !V.allFinite())
		throw std::runtime_error("No finite isosurface; check the full-cell skeleton, radius and mesh resolution");
	Eigen::Vector3d minV = V.colwise().minCoeff();
	Eigen::Vector3d maxV = V.colwise().maxCoeff();

	Eigen::Vector3d move = minV + (maxV - minV) / 2.0;
	for (int j = 0; j < V.rows(); j++)
	{
		V(j, 0) -= move(0);
		V(j, 1) -= move(1);
		V(j, 2) -= move(2);
	}

	//缩放为[-1,-1,-1]到[1,1,1]的边长为2的立方体
	Eigen::Vector3d factor = (maxV - minV) / 2.0;
	double xfac = (1.0 / factor(0));
	double yfac = (1.0 / factor(1));
	double zfac = (1.0 / factor(2));
	for (int j = 0; j < V.rows(); j++)
	{

		V(j, 0) *= xfac;
		V(j, 1) *= yfac;
		V(j, 2) *= zfac;
	}
}

bool Truss::setSolidVolumeFractionModifiedETH(double vof)
{
	//计算当前体积分数
	genSolidTrussVoxelModifiedETH();
	double vof_1 = (double)voxel.sum() / (double)voxel.rows();
	std::cout << "初始体积分数为" << vof_1 << ",半径为" << r_outer << std::endl;

	//将d_outer加一点点,然后计算斜率
	setROuter(r_outer + 0.5, -1);
	genSolidTrussVoxelModifiedETH();
	double vof_2 = (double)voxel.sum() / (double)voxel.rows();
	std::cout << "初始体积分数为" << vof_2 << ",半径为" << r_outer << std::endl;

	double k = ((vof_2 - vof_1) * 1.0) / 0.5;//斜率
	if (!std::isfinite(k)) return false;
	double ini_r = std::abs(k) < 1e-14 ? r_outer : ((vof - vof_2) / k + r_outer);//粗略计算目标相对密度的半径
	if (ini_r < 0)
		ini_r = 1.0;		//从最小开始

	//迭代求解精确的半径,确保体积分数在+-0.5%
	setROuter(ini_r, -1);
	genSolidTrussVoxelModifiedETH();
	double vof_current = (double)voxel.sum() / (double)voxel.rows();
	std::cout << "迭代开始时体积分数为" << vof_current << ",半径为" << ini_r << std::endl;

	double factor = 1.0;	//步进长度因子
	int fit_iterations = 0;
	while (std::abs(vof - vof_current) > 0.005)
	{
		if (factor < 1e-10 || ++fit_iterations > 1000) return false;
		if (vof_current < vof)
		{
			//一点点增加
			while (vof - vof_current > 0.005)
			{
				if (++fit_iterations > 1000) return false;
				ini_r += 0.5 * factor;
				setROuter(ini_r, -1);
				genSolidTrussVoxelModifiedETH();
				vof_current = (double)voxel.sum() / (double)voxel.rows();
			}
		}
		else if (vof_current > vof)
		{
			//一点点减少
			while (vof - vof_current < -0.005)
			{
				if (++fit_iterations > 1000) return false;
				ini_r -= 0.5 * factor;
				setROuter(ini_r, -1);
				genSolidTrussVoxelModifiedETH();
				vof_current = (double)voxel.sum() / (double)voxel.rows();
			}
		}
		else {}
		factor = factor * 0.1;
		std::cout << "迭代时体积分数为" << vof_current << ",半径为" << ini_r << std::endl;
	}
	std::cout << "最终体积分数为" << vof_current << ",半径为" << ini_r << std::endl;
	if (ini_r < 0.8)
		return false;
	return true;
}

void Truss::genSolidTrussVoxelModifiedETH()
{
	double cell_length_half = cell_length / 2;	//模型大小的一半,与一些数值参数有关

	//重建范围
	Eigen::RowVector3d VminCut = { -1 * cell_length_half, -1 * cell_length_half, -1 * cell_length_half };
	Eigen::RowVector3d VmaxCut = { cell_length_half, cell_length_half, cell_length_half };

	Eigen::RowVector3i res = (resolution * ((VmaxCut - VminCut) / (VmaxCut - VminCut).maxCoeff())).cast<int>();

	grid.resize(res(0) * res(1) * res(2), 3);

	//构建采样点的坐标
	const auto lerpCut = [&](const int di, const int d)->double
	{return VminCut(d) + (double)di / (double)(res(d) - 1) * (VmaxCut(d) - VminCut(d)); };

#pragma omp parallel for
	for (int zi = 0; zi < res(2); zi++)
	{
		const double zc = lerpCut(zi, 2);
		for (int yi = 0; yi < res(1); yi++)
		{
			const double yc = lerpCut(yi, 1);
			for (int xi = 0; xi < res(0); xi++)
			{
				const double xc = lerpCut(xi, 0);
				grid.row(xi + res(0) * (yi + res(1) * zi)) = Eigen::RowVector3d(xc, yc, zc);
			}
		}
	}

	//std::vector<Edge> edge_list = base_edges;
	std::vector<Edge> edge_list = edges;
	std::vector<Joint> joint_list = nodes;
	//double smf = smooth_factor;

	//计算距离场的公式
	//看PPL的计算
	//将采样点xyz映射到子单元1/8或者1/48中,然后计算距离,只需要关照子单元中的几条杆
	const auto& function = [edge_list, joint_list, cell_length_half](double x, double y, double z)->double {
		//超出边界,赋一个比较大的正值
		if (abs(x) > cell_length_half || abs(y) > cell_length_half || abs(z) > cell_length_half)
		{
			return 100;
		}
		//将采样点映射到子单元
		Eigen::Vector3d p = { x,y,z };
		getETHSymmetryPoint(p);

		std::vector<double> result;	//存放每个joint的KS函数值,到最后再套一层总的KS函数

		//std::cout << joint_list.size() << std::endl;
		for (int i = 0; i < joint_list.size(); i++)
		{
			double result_joint = 0;		//当前joint的KS函数值
			double smf = joint_list[i].smooth_factor;
			for (int j = 0; j < joint_list[i].associated_edges.size(); j++)
			{
				int index = joint_list[i].associated_edges[j];	//相邻边在edge_list中的索引
				Eigen::Vector3d v1 = edge_list[index].v.row(0), v2 = edge_list[index].v.row(1);
				double dis = pointToSegDist(v1, v2, p) - edge_list[index].r_outer;
				result_joint += exp(-1 / smf * dis);	//smf是公式中的1/s_j
			}
			result_joint = -smf * log(result_joint);

			//判断该joint的n条边是否在空间同一侧,如果是则继续,否则正常计算KS函数值
			if (joint_list[i].is_same_side)
			{
				//判断点p是否在内侧
				Eigen::Vector3d pp = p - joint_list[i].center;
				int n = joint_list[i].associated_edges.size();
				Eigen::MatrixXd E(3, n);

				for (int j = 0; j < n; j++)
				{
					int index = joint_list[i].associated_edges[j];	//相邻边在edge_list中的索引
					Eigen::Vector3d v1 = edge_list[index].v.row(0), v2 = edge_list[index].v.row(1);
					if (isEquation(joint_list[i].center, v1))
						E.col(j) = v2 - joint_list[i].center;
					else if (isEquation(joint_list[i].center, v2))
						E.col(j) = v1 - joint_list[i].center;
					else
						std::cout << "Error in getTopology() calculate is_same_side" << std::endl;
				}
				bool isIn = checkNonNegativeCombination(E, pp);
				//如果在外侧,需要切割
				if (!isIn)
				{
					//计算点p对应哪个边界,即计算到所有扇形面距离减去r,取min
					double minDis = 1e9;
					for (int j = 0; j < n - 1; j++)
					{
						Eigen::Vector3d ej;
						int index_j = joint_list[i].associated_edges[j];	//相邻边在edge_list中的索引
						Eigen::Vector3d v1 = edge_list[index_j].v.row(0), v2 = edge_list[index_j].v.row(1);
						if (isEquation(joint_list[i].center, v1))
							ej = v2;
						else if (isEquation(joint_list[i].center, v2))
							ej = v1;
						else
							std::cout << "Error in getTopology() calculate is_same_side" << std::endl;
						for (int k = j + 1; k < n; k++)
						{
							Eigen::Vector3d ek;
							int index_k = joint_list[i].associated_edges[k];	//相邻边在edge_list中的索引
							Eigen::Vector3d v1 = edge_list[index_k].v.row(0), v2 = edge_list[index_k].v.row(1);
							if (isEquation(joint_list[i].center, v1))
								ek = v2;
							else if (isEquation(joint_list[i].center, v2))
								ek = v1;
							else
								std::cout << "Error in getTopology() calculate is_same_side" << std::endl;

							double dis = pointToTriangleDistance(p, joint_list[i].center, ej, ek);
							dis = dis - std::max(edge_list[index_j].r_outer, edge_list[index_k].r_outer);
							if (dis < minDis)
								minDis = dis;

							/*std::pair<double, double> dis_pair = pointToPartialPlaneDistance(p, joint_list[i].center, ej, ek);
							double dis1 = dis_pair.first - std::max(edge_list[index_j].r_outer, edge_list[index_k].r_outer);
							double dis2 = dis_pair.second - std::max(edge_list[index_j].r_outer, edge_list[index_k].r_outer);
							if (dis1 < minDis1)
								minDis1 = dis1;
							if (dis2 < minDis2)
								minDis2 = dis2;*/
						}
					}

					if (minDis >= 0)
					{
						result_joint = std::max(result_joint, minDis);
					}
					else if (result_joint < 0)
						result_joint = result_joint * (-1) * tanh(minDis);

				}
			}
			result.push_back(result_joint);
		}

		double sum_KS = 0;	//总的KS函数值
		double sum_smf = 0.1;
		for (int i = 0; i < result.size(); i++)
		{
			sum_KS += exp(-1 / sum_smf * result[i]);
		}
		sum_KS = -sum_smf * log(sum_KS);
		return sum_KS;
	};
	voxel.resize(grid.rows());
	voxel.setZero();
#pragma omp parallel for
	for (int i = 0; i < grid.rows(); i++) {
		double indicator = function(grid(i, 0), grid(i, 1), grid(i, 2));
		if (indicator <= 0)
		{
			voxel(i) = 1;
		}
		else
			voxel(i) = 0;
	}
}

void Truss::genSolidTrussMeshModifiedETH()
{
	double cell_length_half = cell_length / 2;	//模型大小的一半,与一些数值参数有关

	//重建范围
	Eigen::RowVector3d Vmin = { -1 * cell_length_half - 0.1, -1 * cell_length_half - 0.1, -1 * cell_length_half - 0.1 };
	Eigen::RowVector3d Vmax = { cell_length_half + 0.1, cell_length_half + 0.1, cell_length_half + 0.1 };

	Eigen::RowVector3i res_mesh = (resolution_mesh * ((Vmax - Vmin) / (Vmax - Vmin).maxCoeff())).cast<int>();

	Eigen::MatrixXd GV(res_mesh(0) * res_mesh(1) * res_mesh(2), 3);

	//构建采样点的坐标
	const auto lerp = [&](const int di, const int d)->double
	{return Vmin(d) + (double)di / (double)(res_mesh(d) - 1) * (Vmax(d) - Vmin(d)); };

#pragma omp parallel for
	for (int zi = 0; zi < res_mesh(2); zi++)
	{
		const double z = lerp(zi, 2);
		for (int yi = 0; yi < res_mesh(1); yi++)
		{
			const double y = lerp(yi, 1);
			for (int xi = 0; xi < res_mesh(0); xi++)
			{
				const double x = lerp(xi, 0);
				GV.row(xi + res_mesh(0) * (yi + res_mesh(1) * zi)) = Eigen::RowVector3d(x, y, z);
			}
		}
	}

	//std::vector<Edge> edge_list = base_edges;
	std::vector<Edge> edge_list = edges;
	std::vector<Joint> joint_list = nodes;
	//double smf = smooth_factor;

	//计算距离场的公式
	//看PPL的计算
	//将采样点xyz映射到子单元1/8或者1/48中,然后计算距离,只需要关照子单元中的几条杆
	const auto& function = [edge_list, joint_list, cell_length_half](double x, double y, double z)->double {
		//超出边界,赋一个比较大的正值
		if (abs(x) > cell_length_half || abs(y) > cell_length_half || abs(z) > cell_length_half)
		{
			return 100;
		}
		//将采样点映射到子单元
		Eigen::Vector3d p = { x,y,z };
		getETHSymmetryPoint(p);

		std::vector<double> result;	//存放每个joint的KS函数值,到最后再套一层总的KS函数

		//std::cout << joint_list.size() << std::endl;
		for (int i = 0; i < joint_list.size(); i++)
		{
			double result_joint = 0;		//当前joint的KS函数值
			double smf = joint_list[i].smooth_factor;
			for (int j = 0; j < joint_list[i].associated_edges.size(); j++)
			{
				int index = joint_list[i].associated_edges[j];	//相邻边在edge_list中的索引
				Eigen::Vector3d v1 = edge_list[index].v.row(0), v2 = edge_list[index].v.row(1);
				double dis = pointToSegDist(v1, v2, p) - edge_list[index].r_outer;
				result_joint += exp(-1 / smf * dis);	//smf是公式中的1/s_j
			}
			result_joint = -smf * log(result_joint);

			//判断该joint的n条边是否在空间同一侧,如果是则继续,否则正常计算KS函数值
			if (joint_list[i].is_same_side)
			{
				//判断点p是否在内侧
				Eigen::Vector3d pp = p - joint_list[i].center;
				int n = joint_list[i].associated_edges.size();
				Eigen::MatrixXd E(3, n);

				for (int j = 0; j < n; j++)
				{
					int index = joint_list[i].associated_edges[j];	//相邻边在edge_list中的索引
					Eigen::Vector3d v1 = edge_list[index].v.row(0), v2 = edge_list[index].v.row(1);
					if (isEquation(joint_list[i].center, v1))
						E.col(j) = v2 - joint_list[i].center;
					else if (isEquation(joint_list[i].center, v2))
						E.col(j) = v1 - joint_list[i].center;
					else
						std::cout << "Error in getTopology() calculate is_same_side" << std::endl;
				}
				bool isIn = checkNonNegativeCombination(E, pp);
				//如果在外侧,需要切割
				if (!isIn)
				{
					//计算点p对应哪个边界,即计算到所有扇形面距离减去r,取min
					double minDis = 1e9;
					for (int j = 0; j < n - 1; j++)
					{
						Eigen::Vector3d ej;
						int index_j = joint_list[i].associated_edges[j];	//相邻边在edge_list中的索引
						Eigen::Vector3d v1 = edge_list[index_j].v.row(0), v2 = edge_list[index_j].v.row(1);
						if (isEquation(joint_list[i].center, v1))
							ej = v2;
						else if (isEquation(joint_list[i].center, v2))
							ej = v1;
						else
							std::cout << "Error in getTopology() calculate is_same_side" << std::endl;
						for (int k = j + 1; k < n; k++)
						{
							Eigen::Vector3d ek;
							int index_k = joint_list[i].associated_edges[k];	//相邻边在edge_list中的索引
							Eigen::Vector3d v1 = edge_list[index_k].v.row(0), v2 = edge_list[index_k].v.row(1);
							if (isEquation(joint_list[i].center, v1))
								ek = v2;
							else if (isEquation(joint_list[i].center, v2))
								ek = v1;
							else
								std::cout << "Error in getTopology() calculate is_same_side" << std::endl;

							double dis = pointToTriangleDistance(p, joint_list[i].center, ej, ek);
							dis = dis - std::max(edge_list[index_j].r_outer, edge_list[index_k].r_outer);
							if (dis < minDis)
								minDis = dis;

							/*std::pair<double, double> dis_pair = pointToPartialPlaneDistance(p, joint_list[i].center, ej, ek);
							double dis1 = dis_pair.first - std::max(edge_list[index_j].r_outer, edge_list[index_k].r_outer);
							double dis2 = dis_pair.second - std::max(edge_list[index_j].r_outer, edge_list[index_k].r_outer);
							if (dis1 < minDis1)
								minDis1 = dis1;
							if (dis2 < minDis2)
								minDis2 = dis2;*/
						}
					}

					if (minDis >= 0)
					{
						result_joint = std::max(result_joint, minDis);
					}
					else if (result_joint < 0)
						result_joint = result_joint * (-1) * tanh(minDis);

				}
			}
			result.push_back(result_joint);
		}

		double sum_KS = 0;	//总的KS函数值
		double sum_smf = 0.1;
		for (int i = 0; i < result.size(); i++)
		{
			sum_KS += exp(-1 / sum_smf * result[i]);
		}
		sum_KS = -sum_smf * log(sum_KS);
		return sum_KS;
	};

	Eigen::VectorXd B(GV.rows());
#pragma omp parallel for
	for (int i = 0; i < GV.rows(); i++) {
		double indicator = function(GV(i, 0), GV(i, 1), GV(i, 2));
		B(i) = indicator;
	}
	igl::copyleft::marching_cubes(B, GV, res_mesh(0), res_mesh(1), res_mesh(2), V, F);

	//平移中心到原点
	if (V.rows() == 0 || !V.allFinite())
		throw std::runtime_error("No finite isosurface; check the full-cell skeleton, radius and mesh resolution");
	Eigen::Vector3d minV = V.colwise().minCoeff();
	Eigen::Vector3d maxV = V.colwise().maxCoeff();

	Eigen::Vector3d move = minV + (maxV - minV) / 2.0;
	for (int j = 0; j < V.rows(); j++)
	{
		V(j, 0) -= move(0);
		V(j, 1) -= move(1);
		V(j, 2) -= move(2);
	}

	//缩放为[-1,-1,-1]到[1,1,1]的边长为2的立方体
	Eigen::Vector3d factor = (maxV - minV) / 2.0;
	double xfac = (1.0 / factor(0));
	double yfac = (1.0 / factor(1));
	double zfac = (1.0 / factor(2));
	for (int j = 0; j < V.rows(); j++)
	{

		V(j, 0) *= xfac;
		V(j, 1) *= yfac;
		V(j, 2) *= zfac;
	}
}

void Truss::genSolidTargetSdfModified(Eigen::MatrixXd& tgrid)
{
	double cell_length_half = cell_length / 2;	//模型大小的一半,与一些数值参数有关

	//重建范围
	Eigen::RowVector3d VminCut = { -1 * cell_length_half, -1 * cell_length_half, -1 * cell_length_half };
	Eigen::RowVector3d VmaxCut = { cell_length_half, cell_length_half, cell_length_half };

	Eigen::RowVector3i res = (resolution * ((VmaxCut - VminCut) / (VmaxCut - VminCut).maxCoeff())).cast<int>();

	grid.resize(res(0) * res(1) * res(2), 3);

	//构建采样点的坐标
	const auto lerpCut = [&](const int di, const int d)->double
	{return VminCut(d) + (double)di / (double)(res(d) - 1) * (VmaxCut(d) - VminCut(d)); };

#pragma omp parallel for
	for (int zi = 0; zi < res(2); zi++)
	{
		const double zc = lerpCut(zi, 2);
		for (int yi = 0; yi < res(1); yi++)
		{
			const double yc = lerpCut(yi, 1);
			for (int xi = 0; xi < res(0); xi++)
			{
				const double xc = lerpCut(xi, 0);
				grid.row(xi + res(0) * (yi + res(1) * zi)) = Eigen::RowVector3d(xc, yc, zc);
			}
		}
	}

	//std::vector<Edge> edge_list = base_edges;
	std::vector<Edge> edge_list = edges;
	std::vector<Joint> joint_list = nodes;
	//double smf = smooth_factor;

	//计算距离场的公式
	//看PPL的计算
	//将采样点xyz映射到子单元1/8或者1/48中,然后计算距离,只需要关照子单元中的几条杆
	const auto& function = [edge_list, joint_list, cell_length_half](double x, double y, double z)->double {
		//超出边界,赋一个比较大的正值
		if (abs(x) > cell_length_half || abs(y) > cell_length_half || abs(z) > cell_length_half)
		{
			return 100;
		}
		//将采样点映射到子单元
		Eigen::Vector3d p = { x,y,z };
		getETHSymmetryPoint(p);

		std::vector<double> result;	//存放每个joint的KS函数值,到最后再套一层总的KS函数

		//std::cout << joint_list.size() << std::endl;
		for (int i = 0; i < joint_list.size(); i++)
		{
			double result_joint = 0;		//当前joint的KS函数值
			double smf = joint_list[i].smooth_factor;
			for (int j = 0; j < joint_list[i].associated_edges.size(); j++)
			{
				int index = joint_list[i].associated_edges[j];	//相邻边在edge_list中的索引
				Eigen::Vector3d v1 = edge_list[index].v.row(0), v2 = edge_list[index].v.row(1);
				double dis = pointToSegDist(v1, v2, p) - edge_list[index].r_outer;
				result_joint += exp(-1 / smf * dis);	//smf是公式中的1/s_j
			}
			result_joint = -smf * log(result_joint);

			//判断该joint的n条边是否在空间同一侧,如果是则继续,否则正常计算KS函数值
			if (joint_list[i].is_same_side)
			{
				//判断点p是否在内侧
				Eigen::Vector3d pp = p - joint_list[i].center;
				int n = joint_list[i].associated_edges.size();
				Eigen::MatrixXd E(3, n);

				for (int j = 0; j < n; j++)
				{
					int index = joint_list[i].associated_edges[j];	//相邻边在edge_list中的索引
					Eigen::Vector3d v1 = edge_list[index].v.row(0), v2 = edge_list[index].v.row(1);
					if (isEquation(joint_list[i].center, v1))
						E.col(j) = v2 - joint_list[i].center;
					else if (isEquation(joint_list[i].center, v2))
						E.col(j) = v1 - joint_list[i].center;
					else
						std::cout << "Error in getTopology() calculate is_same_side" << std::endl;
				}
				bool isIn = checkNonNegativeCombination(E, pp);
				//如果在外侧,需要切割
				if (!isIn)
				{
					//计算点p对应哪个边界,即计算到所有扇形面距离减去r,取min
					double minDis = 1e9;
					for (int j = 0; j < n - 1; j++)
					{
						Eigen::Vector3d ej;
						int index_j = joint_list[i].associated_edges[j];	//相邻边在edge_list中的索引
						Eigen::Vector3d v1 = edge_list[index_j].v.row(0), v2 = edge_list[index_j].v.row(1);
						if (isEquation(joint_list[i].center, v1))
							ej = v2;
						else if (isEquation(joint_list[i].center, v2))
							ej = v1;
						else
							std::cout << "Error in getTopology() calculate is_same_side" << std::endl;
						for (int k = j + 1; k < n; k++)
						{
							Eigen::Vector3d ek;
							int index_k = joint_list[i].associated_edges[k];	//相邻边在edge_list中的索引
							Eigen::Vector3d v1 = edge_list[index_k].v.row(0), v2 = edge_list[index_k].v.row(1);
							if (isEquation(joint_list[i].center, v1))
								ek = v2;
							else if (isEquation(joint_list[i].center, v2))
								ek = v1;
							else
								std::cout << "Error in getTopology() calculate is_same_side" << std::endl;

							double dis = pointToTriangleDistance(p, joint_list[i].center, ej, ek);
							dis = dis - std::max(edge_list[index_j].r_outer, edge_list[index_k].r_outer);
							if (dis < minDis)
								minDis = dis;

							/*std::pair<double, double> dis_pair = pointToPartialPlaneDistance(p, joint_list[i].center, ej, ek);
							double dis1 = dis_pair.first - std::max(edge_list[index_j].r_outer, edge_list[index_k].r_outer);
							double dis2 = dis_pair.second - std::max(edge_list[index_j].r_outer, edge_list[index_k].r_outer);
							if (dis1 < minDis1)
								minDis1 = dis1;
							if (dis2 < minDis2)
								minDis2 = dis2;*/
						}
					}

					if (minDis >= 0)
					{
						result_joint = std::max(result_joint, minDis);
					}
					else if (result_joint < 0)
						result_joint = result_joint * (-1) * tanh(minDis);

				}
			}
			result.push_back(result_joint);
		}

		double sum_KS = 0;	//总的KS函数值
		double sum_smf = 0.1;
		for (int i = 0; i < result.size(); i++)
		{
			sum_KS += exp(-1 / sum_smf * result[i]);
		}
		sum_KS = -sum_smf * log(sum_KS);
		return sum_KS;
	};

	outer_sdf.setZero(tgrid.rows());
	//inner_sdf.setZero(tgrid.rows());

#pragma omp parallel for
	for (int i = 0; i < tgrid.rows(); i++)
	{
		double indicator = function(tgrid(i, 0), tgrid(i, 1), tgrid(i, 2));
		outer_sdf(i) = indicator;
		//实心杆不用inner_sdf
	}

	voxel.resize(grid.rows());
	voxel.setZero();
#pragma omp parallel for
	for (int i = 0; i < grid.rows(); i++) {
		double indicator = function(grid(i, 0), grid(i, 1), grid(i, 2));
		if (indicator <= 0)
		{
			voxel(i) = 1;
		}
		else
			voxel(i) = 0;
	}

	tvol = voxel.sum() / (double)voxel.rows();
}


void Truss::getSolidVelocitiesModified(std::vector<Eigen::VectorXd>& velosities)
{
	//目前设置delta=0.7 可以自行调整
	double delta = 0.7 * (grid.row(1) - grid.row(0)).x();
	double delta2 = pow(delta, 2);	//没用上
	delta = delta2;
	//double coef = 0.75 / delta;		//没用上

	Truss PTL2 = *this;
	double perturbation_R = 1.0, perturbation_smf = 0.3;		//步进值,需要自行调整

	//计算当前的体素和距离场
	Eigen::MatrixXd tgrid;	//值为1的体素坐标
	Eigen::VectorXi map;
	getTargetCoors(tgrid, map);
	genSolidTargetSdfModified(tgrid);
	int vsize = tgrid.rows();

	velosities.clear();
	dV.clear();
	dH.clear();

	//计算R改变后的体素和距离场,并计算差值
	for (int index = 0; index < base_edges.size(); index++)
	{
		//int index = 0;
		Eigen::VectorXd vv;
		double current_r_outer = PTL2.getROuter(index);
		PTL2.setROuter(current_r_outer + perturbation_R, index);
		PTL2.genSolidTargetSdfModified(tgrid);
		PTL2.setROuter(current_r_outer, index);	//复原

		vv.resize(vsize);
		vv.setZero();
		Eigen::VectorXd sdf_diff;	//两个距离场的差
		sdf_diff.setZero(vsize);
		Eigen::VectorXd outer_sdf_diff = outer_sdf - PTL2.outer_sdf;
		//Eigen::VectorXd inner_sdf_diff = inner_sdf - s2.inner_sdf;

		Eigen::VectorXi vdiff = PTL2.voxel - voxel;	//为什么这里是PTL2-PTL1? 因为距离场diff取了绝对值

		Eigen::VectorXd diffp(vsize);
		Eigen::VectorXd diffn(vsize);
		Eigen::VectorXi vdiffp(voxel.rows());
		Eigen::VectorXi vdiffn(voxel.rows());
		diffp.setZero();
		diffn.setZero();
		vdiffp.setZero();
		vdiffn.setZero();

#pragma omp parallel for
		for (int j = 0; j < voxel.rows(); j++) {
			if (voxel(j) == 1) {
				//需要是边界上的体素点
				if ((!isInnerVoxel(j))) {
					int k = map(j);	//取出在tgrid中的坐标

					//Dirac function
					//identify outer boundary
					if (abs(outer_sdf(k)) <= delta) {
						sdf_diff(k) = outer_sdf_diff(k);
					}
					//inner boundary
					//else if (abs(inner_sdf(k) <= delta)) {
					//	sdf_diff(k) = inner_sdf_diff(k);
					//}
					else {
						sdf_diff(k) = 0;
					}
					if (sdf_diff(k) > 0) {
						diffp(k) = sdf_diff(k);
					}
					if (sdf_diff(k) < 0) {
						diffn(k) = sdf_diff(k);
					}
				}
			}
			if (vdiff(j) > 0) {
				vdiffp(j) = vdiff(j);
			}
			else if (vdiff(j) < 0) {
				vdiffn(j) = vdiff(j);
			}
		}

		double sum_diffp = diffp.sum();
		double sum_diffn = diffn.sum();
		double sum_vdiffp = vdiffp.sum();
		double sum_vdiffn = vdiffn.sum();

		double factor1, factor2;
		if (sum_diffp == 0) {
			factor1 = 0;
		}
		else {
			factor1 = sum_vdiffp / sum_diffp;
		}
		if (sum_diffn == 0) {

			factor2 = 0;
		}
		else {
			factor2 = sum_vdiffn / sum_diffn;
		}

		vv = factor1 * diffp + factor2 * diffn;
		vv = vv / perturbation_R;

		double dv = (PTL2.tvol - tvol) / perturbation_R;
		velosities.push_back(vv);
		dV.push_back(dv);
	}
	//计算smf改变后的体素和距离场,并计算差值
	for (int index = 0; index < base_joints.size(); index++)
	{
		Eigen::VectorXd vv;
		double current_smooth_factor = PTL2.getSmoothFactor(index);
		PTL2.setSmoothFactor(current_smooth_factor + perturbation_smf, index);
		PTL2.genSolidTargetSdfModified(tgrid);
		PTL2.setSmoothFactor(current_smooth_factor, index);	//复原

		vv.resize(vsize);
		vv.setZero();
		Eigen::VectorXd sdf_diff;	//两个距离场的差
		sdf_diff.setZero(vsize);
		Eigen::VectorXd outer_sdf_diff = outer_sdf - PTL2.outer_sdf;
		//Eigen::VectorXd inner_sdf_diff = inner_sdf - s2.inner_sdf;

		Eigen::VectorXi vdiff = PTL2.voxel - voxel;	//为什么这里是PTL2-PTL1? 因为距离场diff取了绝对值

		Eigen::VectorXd diffp(vsize);
		Eigen::VectorXd diffn(vsize);
		Eigen::VectorXi vdiffp(voxel.rows());
		Eigen::VectorXi vdiffn(voxel.rows());
		diffp.setZero();
		diffn.setZero();
		vdiffp.setZero();
		vdiffn.setZero();

#pragma omp parallel for
		for (int j = 0; j < voxel.rows(); j++) {
			if (voxel(j) == 1) {
				//需要是边界上的体素点
				if ((!isInnerVoxel(j))) {
					int k = map(j);	//取出在tgrid中的坐标

					//Dirac function
					//identify outer boundary
					if (abs(outer_sdf(k)) <= delta) {
						sdf_diff(k) = outer_sdf_diff(k);
					}
					//inner boundary
					//else if (abs(inner_sdf(k) <= delta)) {
					//	sdf_diff(k) = inner_sdf_diff(k);
					//}
					else {
						sdf_diff(k) = 0;
					}
					if (sdf_diff(k) > 0) {
						diffp(k) = sdf_diff(k);
					}
					if (sdf_diff(k) < 0) {
						diffn(k) = sdf_diff(k);
					}
				}
			}
			if (vdiff(j) > 0) {
				vdiffp(j) = vdiff(j);
			}
			else if (vdiff(j) < 0) {
				vdiffn(j) = vdiff(j);
			}
		}

		double sum_diffp = diffp.sum();
		double sum_diffn = diffn.sum();
		double sum_vdiffp = vdiffp.sum();
		double sum_vdiffn = vdiffn.sum();

		double factor1, factor2;
		if (sum_diffp == 0) {
			factor1 = 0;
		}
		else {
			factor1 = sum_vdiffp / sum_diffp;
		}
		if (sum_diffn == 0) {

			factor2 = 0;
		}
		else {
			factor2 = sum_vdiffn / sum_diffn;
		}

		vv = factor1 * diffp + factor2 * diffn;
		vv = vv / perturbation_smf;

		double dv = (PTL2.tvol - tvol) / perturbation_smf;
		velosities.push_back(vv);
		dV.push_back(dv);
	}
}

void Truss::genTilingTrussMeshModified(int size_x, int size_y, int size_z)
{
	//应注意重建范围 [0,0,0]->[cell_length*size_x,cell_length*size_y,cell_length*size_z] 再加入边界一小圈
	Eigen::RowVector3d Vmin = { 0.0, 0.0, 0.0 };
	Eigen::RowVector3d Vmax = { cell_length * size_x, cell_length * size_y, cell_length * size_z };

	//这里减掉一些是因为每两个格子会重复用一个采样点
	Eigen::RowVector3i res_mesh = { resolution_mesh * size_x - (size_x - 1),resolution_mesh * size_y - (size_y - 1),resolution_mesh * size_z - (size_z - 1) };

	Eigen::MatrixXd GV((res_mesh(0) + 2) * (res_mesh(1) + 2) * (res_mesh(2) + 2), 3);	//算上边界

	//构建采样点的坐标
	const auto lerp = [&](const int di, const int d)->double
	{return Vmin(d) + (double)di / (double)(res_mesh(d) - 1) * (Vmax(d) - Vmin(d)); };

	double boundary_length = 0.1;

#pragma omp parallel for
	for (int zi = 0; zi < res_mesh(2) + 2; zi++)
	{
		double z = 0;
		if (zi == 0)
			z = Vmin(2) - boundary_length;
		else if (zi == res_mesh(2) + 1)
			z = Vmax(2) + boundary_length;
		else
			z = lerp(zi, 2);
		for (int yi = 0; yi < res_mesh(1) + 2; yi++)
		{
			double y = 0;
			if (yi == 0)
				y = Vmin(1) - boundary_length;
			else if (yi == res_mesh(1) + 1)
				y = Vmax(1) + boundary_length;
			else
				y = lerp(yi, 1);
			for (int xi = 0; xi < res_mesh(0) + 2; xi++)
			{
				double x = 0;
				if (xi == 0)
					x = Vmin(0) - boundary_length;
				else if (xi == res_mesh(0) + 1)
					x = Vmax(0) + boundary_length;
				else
					x = lerp(xi, 0);
				GV.row(xi + (res_mesh(0) + 2) * (yi + (res_mesh(1) + 2) * zi)) = Eigen::RowVector3d(x, y, z);
			}
		}
	}

	//std::vector<Edge> edge_list = base_edges;
	std::vector<Edge> edge_list = edges;
	std::vector<Joint> joint_list = nodes;
	double cell_len = cell_length;
	//double smf = smooth_factor;

	//计算距离场的公式
	//看PPL的计算
	//将采样点xyz映射到子单元1/8或者1/48中,然后计算距离,只需要关照子单元中的几条杆
	const auto& function = [edge_list, joint_list, Vmin, Vmax, cell_len](double x, double y, double z)->double {
		//超出边界,赋一个比较大的正值
		if (x<Vmin(0) || x>Vmax(0) || y<Vmin(1) || y>Vmax(1) || z<Vmin(2) || z>Vmax(2))
		{
			return 100;
		}
		//将采样点映射到子单元,这里注意先取余到第一个单元[0,0,0]->[cell_length,cell_length,cell_length],再平移到以0为原点
		while (x > cell_len)
			x = x - cell_len;
		while (y > cell_len)
			y = y - cell_len;
		while (z > cell_len)
			z = z - cell_len;
		Eigen::Vector3d p = { x - cell_len / 2.0,y - cell_len / 2.0,z - cell_len / 2.0 };
		getPanettaSymmetryPoint(p);

		std::vector<double> result;	//存放每个joint的KS函数值,到最后再套一层总的KS函数

		for (int i = 0; i < joint_list.size(); i++)
		{
			double result_joint = 0;		//当前joint的KS函数值
			double smf = joint_list[i].smooth_factor;
			for (int j = 0; j < joint_list[i].associated_edges.size(); j++)
			{
				int index = joint_list[i].associated_edges[j];	//相邻边在edge_list中的索引
				Eigen::Vector3d v1 = edge_list[index].v.row(0), v2 = edge_list[index].v.row(1);
				double dis = pointToSegDist(v1, v2, p) - edge_list[index].r_outer;
				result_joint += exp(-1 / smf * dis);	//smf是公式中的1/s_j
			}
			result_joint = -smf * log(result_joint);

			//判断该joint的n条边是否在空间同一侧,如果是则继续,否则正常计算KS函数值
			if (joint_list[i].is_same_side)
			{
				//判断点p是否在内侧
				Eigen::Vector3d pp = p - joint_list[i].center;
				int n = joint_list[i].associated_edges.size();
				Eigen::MatrixXd E(3, n);

				for (int j = 0; j < n; j++)
				{
					int index = joint_list[i].associated_edges[j];	//相邻边在edge_list中的索引
					Eigen::Vector3d v1 = edge_list[index].v.row(0), v2 = edge_list[index].v.row(1);
					if (isEquation(joint_list[i].center, v1))
						E.col(j) = v2 - joint_list[i].center;
					else if (isEquation(joint_list[i].center, v2))
						E.col(j) = v1 - joint_list[i].center;
					else
						std::cout << "Error in getTopology() calculate is_same_side" << std::endl;
				}
				bool isIn = checkNonNegativeCombination(E, pp);
				//如果在外侧,需要切割
				if (!isIn)
				{
					//计算点p对应哪个边界,即计算到所有扇形面距离减去r,取min
					double minDis = 1e9;
					for (int j = 0; j < n - 1; j++)
					{
						Eigen::Vector3d ej;
						int index_j = joint_list[i].associated_edges[j];	//相邻边在edge_list中的索引
						Eigen::Vector3d v1 = edge_list[index_j].v.row(0), v2 = edge_list[index_j].v.row(1);
						if (isEquation(joint_list[i].center, v1))
							ej = v2;
						else if (isEquation(joint_list[i].center, v2))
							ej = v1;
						else
							std::cout << "Error in getTopology() calculate is_same_side" << std::endl;
						for (int k = j + 1; k < n; k++)
						{
							Eigen::Vector3d ek;
							int index_k = joint_list[i].associated_edges[k];	//相邻边在edge_list中的索引
							Eigen::Vector3d v1 = edge_list[index_k].v.row(0), v2 = edge_list[index_k].v.row(1);
							if (isEquation(joint_list[i].center, v1))
								ek = v2;
							else if (isEquation(joint_list[i].center, v2))
								ek = v1;
							else
								std::cout << "Error in getTopology() calculate is_same_side" << std::endl;

							double dis = pointToTriangleDistance(p, joint_list[i].center, ej, ek);
							dis = dis - std::max(edge_list[index_j].r_outer, edge_list[index_k].r_outer);
							if (dis < minDis)
								minDis = dis;

							/*std::pair<double, double> dis_pair = pointToPartialPlaneDistance(p, joint_list[i].center, ej, ek);
							double dis1 = dis_pair.first - std::max(edge_list[index_j].r_outer, edge_list[index_k].r_outer);
							double dis2 = dis_pair.second - std::max(edge_list[index_j].r_outer, edge_list[index_k].r_outer);
							if (dis1 < minDis1)
								minDis1 = dis1;
							if (dis2 < minDis2)
								minDis2 = dis2;*/
						}
					}

					if (minDis >= 0)
					{
						result_joint = std::max(result_joint, minDis);
					}
					else if (result_joint < 0)
						result_joint = result_joint * (-1) * tanh(minDis);

				}
			}
			result.push_back(result_joint);
		}

		double sum_KS = 0;	//总的KS函数值
		double sum_smf = 0.1;
		for (int i = 0; i < result.size(); i++)
		{
			sum_KS += exp(-1 / sum_smf * result[i]);
		}
		sum_KS = -sum_smf * log(sum_KS);
		return sum_KS;
	};


	Eigen::VectorXd B(GV.rows());
#pragma omp parallel for
	for (int i = 0; i < GV.rows(); i++) {
		double indicator = function(GV(i, 0), GV(i, 1), GV(i, 2));
		B(i) = indicator;
	}
	igl::copyleft::marching_cubes(B, GV, (res_mesh(0) + 2), (res_mesh(1) + 2), (res_mesh(2) + 2), V, F);

	//平移min到原点
	if (V.rows() == 0 || !V.allFinite())
		throw std::runtime_error("No finite isosurface; check the full-cell skeleton, radius and mesh resolution");
	Eigen::Vector3d minV = V.colwise().minCoeff();
	Eigen::Vector3d maxV = V.colwise().maxCoeff();

	//Eigen::Vector3d move = minV + (maxV - minV) / 2.0;
	Eigen::Vector3d move = minV;
	for (int j = 0; j < V.rows(); j++)
	{
		V(j, 0) -= move(0);
		V(j, 1) -= move(1);
		V(j, 2) -= move(2);
	}

	//缩放为[0,0,0]开始的边长为1的立方体
	Eigen::Vector3d factor = maxV - minV;
	double xfac = ((double)size_x / factor(0));
	double yfac = ((double)size_y / factor(1));
	double zfac = ((double)size_z / factor(2));
	for (int j = 0; j < V.rows(); j++)
	{

		V(j, 0) *= xfac;
		V(j, 1) *= yfac;
		V(j, 2) *= zfac;
	}
}

void Truss::genSolidTargetSdfModifiedETH(Eigen::MatrixXd& tgrid)
{
	double cell_length_half = cell_length / 2;	//模型大小的一半,与一些数值参数有关

	//重建范围
	Eigen::RowVector3d VminCut = { -1 * cell_length_half, -1 * cell_length_half, -1 * cell_length_half };
	Eigen::RowVector3d VmaxCut = { cell_length_half, cell_length_half, cell_length_half };

	Eigen::RowVector3i res = (resolution * ((VmaxCut - VminCut) / (VmaxCut - VminCut).maxCoeff())).cast<int>();

	grid.resize(res(0) * res(1) * res(2), 3);

	//构建采样点的坐标
	const auto lerpCut = [&](const int di, const int d)->double
	{return VminCut(d) + (double)di / (double)(res(d) - 1) * (VmaxCut(d) - VminCut(d)); };

#pragma omp parallel for
	for (int zi = 0; zi < res(2); zi++)
	{
		const double zc = lerpCut(zi, 2);
		for (int yi = 0; yi < res(1); yi++)
		{
			const double yc = lerpCut(yi, 1);
			for (int xi = 0; xi < res(0); xi++)
			{
				const double xc = lerpCut(xi, 0);
				grid.row(xi + res(0) * (yi + res(1) * zi)) = Eigen::RowVector3d(xc, yc, zc);
			}
		}
	}

	//std::vector<Edge> edge_list = base_edges;
	std::vector<Edge> edge_list = edges;
	std::vector<Joint> joint_list = nodes;
	//double smf = smooth_factor;

	//计算距离场的公式
	//看PPL的计算
	//将采样点xyz映射到子单元1/8或者1/48中,然后计算距离,只需要关照子单元中的几条杆
	const auto& function = [edge_list, joint_list, cell_length_half](double x, double y, double z)->double {
		//超出边界,赋一个比较大的正值
		if (abs(x) > cell_length_half || abs(y) > cell_length_half || abs(z) > cell_length_half)
		{
			return 100;
		}
		//将采样点映射到子单元
		Eigen::Vector3d p = { x,y,z };
		getETHSymmetryPoint(p);

		std::vector<double> result;	//存放每个joint的KS函数值,到最后再套一层总的KS函数

		//std::cout << joint_list.size() << std::endl;
		for (int i = 0; i < joint_list.size(); i++)
		{
			double result_joint = 0;		//当前joint的KS函数值
			double smf = joint_list[i].smooth_factor;
			for (int j = 0; j < joint_list[i].associated_edges.size(); j++)
			{
				int index = joint_list[i].associated_edges[j];	//相邻边在edge_list中的索引
				Eigen::Vector3d v1 = edge_list[index].v.row(0), v2 = edge_list[index].v.row(1);
				double dis = pointToSegDist(v1, v2, p) - edge_list[index].r_outer;
				result_joint += exp(-1 / smf * dis);	//smf是公式中的1/s_j
			}
			result_joint = -smf * log(result_joint);

			//判断该joint的n条边是否在空间同一侧,如果是则继续,否则正常计算KS函数值
			if (joint_list[i].is_same_side)
			{
				//判断点p是否在内侧
				Eigen::Vector3d pp = p - joint_list[i].center;
				int n = joint_list[i].associated_edges.size();
				Eigen::MatrixXd E(3, n);

				for (int j = 0; j < n; j++)
				{
					int index = joint_list[i].associated_edges[j];	//相邻边在edge_list中的索引
					Eigen::Vector3d v1 = edge_list[index].v.row(0), v2 = edge_list[index].v.row(1);
					if (isEquation(joint_list[i].center, v1))
						E.col(j) = v2 - joint_list[i].center;
					else if (isEquation(joint_list[i].center, v2))
						E.col(j) = v1 - joint_list[i].center;
					else
						std::cout << "Error in getTopology() calculate is_same_side" << std::endl;
				}
				bool isIn = checkNonNegativeCombination(E, pp);
				//如果在外侧,需要切割
				if (!isIn)
				{
					//计算点p对应哪个边界,即计算到所有扇形面距离减去r,取min
					double minDis = 1e9;
					for (int j = 0; j < n - 1; j++)
					{
						Eigen::Vector3d ej;
						int index_j = joint_list[i].associated_edges[j];	//相邻边在edge_list中的索引
						Eigen::Vector3d v1 = edge_list[index_j].v.row(0), v2 = edge_list[index_j].v.row(1);
						if (isEquation(joint_list[i].center, v1))
							ej = v2;
						else if (isEquation(joint_list[i].center, v2))
							ej = v1;
						else
							std::cout << "Error in getTopology() calculate is_same_side" << std::endl;
						for (int k = j + 1; k < n; k++)
						{
							Eigen::Vector3d ek;
							int index_k = joint_list[i].associated_edges[k];	//相邻边在edge_list中的索引
							Eigen::Vector3d v1 = edge_list[index_k].v.row(0), v2 = edge_list[index_k].v.row(1);
							if (isEquation(joint_list[i].center, v1))
								ek = v2;
							else if (isEquation(joint_list[i].center, v2))
								ek = v1;
							else
								std::cout << "Error in getTopology() calculate is_same_side" << std::endl;

							double dis = pointToTriangleDistance(p, joint_list[i].center, ej, ek);
							dis = dis - std::max(edge_list[index_j].r_outer, edge_list[index_k].r_outer);
							if (dis < minDis)
								minDis = dis;

							/*std::pair<double, double> dis_pair = pointToPartialPlaneDistance(p, joint_list[i].center, ej, ek);
							double dis1 = dis_pair.first - std::max(edge_list[index_j].r_outer, edge_list[index_k].r_outer);
							double dis2 = dis_pair.second - std::max(edge_list[index_j].r_outer, edge_list[index_k].r_outer);
							if (dis1 < minDis1)
								minDis1 = dis1;
							if (dis2 < minDis2)
								minDis2 = dis2;*/
						}
					}

					if (minDis >= 0)
					{
						result_joint = std::max(result_joint, minDis);
					}
					else if (result_joint < 0)
						result_joint = result_joint * (-1) * tanh(minDis);

				}
			}
			result.push_back(result_joint);
		}

		double sum_KS = 0;	//总的KS函数值
		double sum_smf = 0.1;
		for (int i = 0; i < result.size(); i++)
		{
			sum_KS += exp(-1 / sum_smf * result[i]);
		}
		sum_KS = -sum_smf * log(sum_KS);
		return sum_KS;
	};

	outer_sdf.setZero(tgrid.rows());
	//inner_sdf.setZero(tgrid.rows());

#pragma omp parallel for
	for (int i = 0; i < tgrid.rows(); i++)
	{
		double indicator = function(tgrid(i, 0), tgrid(i, 1), tgrid(i, 2));
		outer_sdf(i) = indicator;
		//实心杆不用inner_sdf
	}

	voxel.resize(grid.rows());
	voxel.setZero();
#pragma omp parallel for
	for (int i = 0; i < grid.rows(); i++) {
		double indicator = function(grid(i, 0), grid(i, 1), grid(i, 2));
		if (indicator <= 0)
		{
			voxel(i) = 1;
		}
		else
			voxel(i) = 0;
	}

	tvol = voxel.sum() / (double)voxel.rows();
}


void Truss::getSolidVelocitiesModifiedETH(std::vector<Eigen::VectorXd>& velosities)
{
	//目前设置delta=0.7 可以自行调整
	double delta = 0.7 * (grid.row(1) - grid.row(0)).x();
	double delta2 = pow(delta, 2);	//没用上
	delta = delta2;
	//double coef = 0.75 / delta;		//没用上

	Truss PTL2 = *this;
	double perturbation_R = 1.0, perturbation_smf = 0.3;		//步进值,需要自行调整

	//计算当前的体素和距离场
	Eigen::MatrixXd tgrid;	//值为1的体素坐标
	Eigen::VectorXi map;
	getTargetCoors(tgrid, map);
	genSolidTargetSdfModifiedETH(tgrid);
	int vsize = tgrid.rows();

	velosities.clear();
	dV.clear();
	dH.clear();

	//计算R改变后的体素和距离场,并计算差值
	for (int index = 0; index < base_edges.size(); index++)
	{
		//int index = 0;
		Eigen::VectorXd vv;
		double current_r_outer = PTL2.getROuter(index);
		PTL2.setROuter(current_r_outer + perturbation_R, index);
		PTL2.genSolidTargetSdfModifiedETH(tgrid);
		PTL2.setROuter(current_r_outer, index);	//复原

		vv.resize(vsize);
		vv.setZero();
		Eigen::VectorXd sdf_diff;	//两个距离场的差
		sdf_diff.setZero(vsize);
		Eigen::VectorXd outer_sdf_diff = outer_sdf - PTL2.outer_sdf;
		//Eigen::VectorXd inner_sdf_diff = inner_sdf - s2.inner_sdf;

		Eigen::VectorXi vdiff = PTL2.voxel - voxel;	//为什么这里是PTL2-PTL1? 因为距离场diff取了绝对值

		Eigen::VectorXd diffp(vsize);
		Eigen::VectorXd diffn(vsize);
		Eigen::VectorXi vdiffp(voxel.rows());
		Eigen::VectorXi vdiffn(voxel.rows());
		diffp.setZero();
		diffn.setZero();
		vdiffp.setZero();
		vdiffn.setZero();

#pragma omp parallel for
		for (int j = 0; j < voxel.rows(); j++) {
			if (voxel(j) == 1) {
				//需要是边界上的体素点
				if ((!isInnerVoxel(j))) {
					int k = map(j);	//取出在tgrid中的坐标

					//Dirac function
					//identify outer boundary
					if (abs(outer_sdf(k)) <= delta) {
						sdf_diff(k) = outer_sdf_diff(k);
					}
					//inner boundary
					//else if (abs(inner_sdf(k) <= delta)) {
					//	sdf_diff(k) = inner_sdf_diff(k);
					//}
					else {
						sdf_diff(k) = 0;
					}
					if (sdf_diff(k) > 0) {
						diffp(k) = sdf_diff(k);
					}
					if (sdf_diff(k) < 0) {
						diffn(k) = sdf_diff(k);
					}
				}
			}
			if (vdiff(j) > 0) {
				vdiffp(j) = vdiff(j);
			}
			else if (vdiff(j) < 0) {
				vdiffn(j) = vdiff(j);
			}
		}

		double sum_diffp = diffp.sum();
		double sum_diffn = diffn.sum();
		double sum_vdiffp = vdiffp.sum();
		double sum_vdiffn = vdiffn.sum();

		double factor1, factor2;
		if (sum_diffp == 0) {
			factor1 = 0;
		}
		else {
			factor1 = sum_vdiffp / sum_diffp;
		}
		if (sum_diffn == 0) {

			factor2 = 0;
		}
		else {
			factor2 = sum_vdiffn / sum_diffn;
		}

		vv = factor1 * diffp + factor2 * diffn;
		vv = vv / perturbation_R;

		double dv = (PTL2.tvol - tvol) / perturbation_R;
		velosities.push_back(vv);
		dV.push_back(dv);
	}
	//计算smf改变后的体素和距离场,并计算差值
	for (int index = 0; index < base_joints.size(); index++)
	{
		Eigen::VectorXd vv;
		double current_smooth_factor = PTL2.getSmoothFactor(index);
		PTL2.setSmoothFactor(current_smooth_factor + perturbation_smf, index);
		PTL2.genSolidTargetSdfModifiedETH(tgrid);
		PTL2.setSmoothFactor(current_smooth_factor, index);	//复原

		vv.resize(vsize);
		vv.setZero();
		Eigen::VectorXd sdf_diff;	//两个距离场的差
		sdf_diff.setZero(vsize);
		Eigen::VectorXd outer_sdf_diff = outer_sdf - PTL2.outer_sdf;
		//Eigen::VectorXd inner_sdf_diff = inner_sdf - s2.inner_sdf;

		Eigen::VectorXi vdiff = PTL2.voxel - voxel;	//为什么这里是PTL2-PTL1? 因为距离场diff取了绝对值

		Eigen::VectorXd diffp(vsize);
		Eigen::VectorXd diffn(vsize);
		Eigen::VectorXi vdiffp(voxel.rows());
		Eigen::VectorXi vdiffn(voxel.rows());
		diffp.setZero();
		diffn.setZero();
		vdiffp.setZero();
		vdiffn.setZero();

#pragma omp parallel for
		for (int j = 0; j < voxel.rows(); j++) {
			if (voxel(j) == 1) {
				//需要是边界上的体素点
				if ((!isInnerVoxel(j))) {
					int k = map(j);	//取出在tgrid中的坐标

					//Dirac function
					//identify outer boundary
					if (abs(outer_sdf(k)) <= delta) {
						sdf_diff(k) = outer_sdf_diff(k);
					}
					//inner boundary
					//else if (abs(inner_sdf(k) <= delta)) {
					//	sdf_diff(k) = inner_sdf_diff(k);
					//}
					else {
						sdf_diff(k) = 0;
					}
					if (sdf_diff(k) > 0) {
						diffp(k) = sdf_diff(k);
					}
					if (sdf_diff(k) < 0) {
						diffn(k) = sdf_diff(k);
					}
				}
			}
			if (vdiff(j) > 0) {
				vdiffp(j) = vdiff(j);
			}
			else if (vdiff(j) < 0) {
				vdiffn(j) = vdiff(j);
			}
		}

		double sum_diffp = diffp.sum();
		double sum_diffn = diffn.sum();
		double sum_vdiffp = vdiffp.sum();
		double sum_vdiffn = vdiffn.sum();

		double factor1, factor2;
		if (sum_diffp == 0) {
			factor1 = 0;
		}
		else {
			factor1 = sum_vdiffp / sum_diffp;
		}
		if (sum_diffn == 0) {

			factor2 = 0;
		}
		else {
			factor2 = sum_vdiffn / sum_diffn;
		}

		vv = factor1 * diffp + factor2 * diffn;
		vv = vv / perturbation_smf;

		double dv = (PTL2.tvol - tvol) / perturbation_smf;
		velosities.push_back(vv);
		dV.push_back(dv);
	}
}

void Truss::genTilingTrussMeshModifiedETH(int size_x, int size_y, int size_z)
{
	//应注意重建范围 [0,0,0]->[cell_length*size_x,cell_length*size_y,cell_length*size_z] 再加入边界一小圈
	Eigen::RowVector3d Vmin = { 0.0, 0.0, 0.0 };
	Eigen::RowVector3d Vmax = { cell_length * size_x, cell_length * size_y, cell_length * size_z };

	//这里减掉一些是因为每两个格子会重复用一个采样点
	Eigen::RowVector3i res_mesh = { resolution_mesh * size_x - (size_x - 1),resolution_mesh * size_y - (size_y - 1),resolution_mesh * size_z - (size_z - 1) };

	Eigen::MatrixXd GV((res_mesh(0) + 2) * (res_mesh(1) + 2) * (res_mesh(2) + 2), 3);	//算上边界

	//构建采样点的坐标
	const auto lerp = [&](const int di, const int d)->double
	{return Vmin(d) + (double)di / (double)(res_mesh(d) - 1) * (Vmax(d) - Vmin(d)); };

	double boundary_length = 0.1;

#pragma omp parallel for
	for (int zi = 0; zi < res_mesh(2) + 2; zi++)
	{
		double z = 0;
		if (zi == 0)
			z = Vmin(2) - boundary_length;
		else if (zi == res_mesh(2) + 1)
			z = Vmax(2) + boundary_length;
		else
			z = lerp(zi, 2);
		for (int yi = 0; yi < res_mesh(1) + 2; yi++)
		{
			double y = 0;
			if (yi == 0)
				y = Vmin(1) - boundary_length;
			else if (yi == res_mesh(1) + 1)
				y = Vmax(1) + boundary_length;
			else
				y = lerp(yi, 1);
			for (int xi = 0; xi < res_mesh(0) + 2; xi++)
			{
				double x = 0;
				if (xi == 0)
					x = Vmin(0) - boundary_length;
				else if (xi == res_mesh(0) + 1)
					x = Vmax(0) + boundary_length;
				else
					x = lerp(xi, 0);
				GV.row(xi + (res_mesh(0) + 2) * (yi + (res_mesh(1) + 2) * zi)) = Eigen::RowVector3d(x, y, z);
			}
		}
	}

	//std::vector<Edge> edge_list = base_edges;
	std::vector<Edge> edge_list = edges;
	std::vector<Joint> joint_list = nodes;
	double cell_len = cell_length;
	//double smf = smooth_factor;

	//计算距离场的公式
	//看PPL的计算
	//将采样点xyz映射到子单元1/8或者1/48中,然后计算距离,只需要关照子单元中的几条杆
	const auto& function = [edge_list, joint_list, Vmin, Vmax, cell_len](double x, double y, double z)->double {
		//超出边界,赋一个比较大的正值
		if (x<Vmin(0) || x>Vmax(0) || y<Vmin(1) || y>Vmax(1) || z<Vmin(2) || z>Vmax(2))
		{
			return 100;
		}
		//将采样点映射到子单元,这里注意先取余到第一个单元[0,0,0]->[cell_length,cell_length,cell_length],再平移到以0为原点
		while (x > cell_len)
			x = x - cell_len;
		while (y > cell_len)
			y = y - cell_len;
		while (z > cell_len)
			z = z - cell_len;
		Eigen::Vector3d p = { x - cell_len / 2.0,y - cell_len / 2.0,z - cell_len / 2.0 };
		getETHSymmetryPoint(p);

		std::vector<double> result;	//存放每个joint的KS函数值,到最后再套一层总的KS函数

		for (int i = 0; i < joint_list.size(); i++)
		{
			double result_joint = 0;		//当前joint的KS函数值
			double smf = joint_list[i].smooth_factor;
			for (int j = 0; j < joint_list[i].associated_edges.size(); j++)
			{
				int index = joint_list[i].associated_edges[j];	//相邻边在edge_list中的索引
				Eigen::Vector3d v1 = edge_list[index].v.row(0), v2 = edge_list[index].v.row(1);
				double dis = pointToSegDist(v1, v2, p) - edge_list[index].r_outer;
				result_joint += exp(-1 / smf * dis);	//smf是公式中的1/s_j
			}
			result_joint = -smf * log(result_joint);

			//判断该joint的n条边是否在空间同一侧,如果是则继续,否则正常计算KS函数值
			if (joint_list[i].is_same_side)
			{
				//判断点p是否在内侧
				Eigen::Vector3d pp = p - joint_list[i].center;
				int n = joint_list[i].associated_edges.size();
				Eigen::MatrixXd E(3, n);

				for (int j = 0; j < n; j++)
				{
					int index = joint_list[i].associated_edges[j];	//相邻边在edge_list中的索引
					Eigen::Vector3d v1 = edge_list[index].v.row(0), v2 = edge_list[index].v.row(1);
					if (isEquation(joint_list[i].center, v1))
						E.col(j) = v2 - joint_list[i].center;
					else if (isEquation(joint_list[i].center, v2))
						E.col(j) = v1 - joint_list[i].center;
					else
						std::cout << "Error in getTopology() calculate is_same_side" << std::endl;
				}
				bool isIn = checkNonNegativeCombination(E, pp);
				//如果在外侧,需要切割
				if (!isIn)
				{
					//计算点p对应哪个边界,即计算到所有扇形面距离减去r,取min
					double minDis = 1e9;
					for (int j = 0; j < n - 1; j++)
					{
						Eigen::Vector3d ej;
						int index_j = joint_list[i].associated_edges[j];	//相邻边在edge_list中的索引
						Eigen::Vector3d v1 = edge_list[index_j].v.row(0), v2 = edge_list[index_j].v.row(1);
						if (isEquation(joint_list[i].center, v1))
							ej = v2;
						else if (isEquation(joint_list[i].center, v2))
							ej = v1;
						else
							std::cout << "Error in getTopology() calculate is_same_side" << std::endl;
						for (int k = j + 1; k < n; k++)
						{
							Eigen::Vector3d ek;
							int index_k = joint_list[i].associated_edges[k];	//相邻边在edge_list中的索引
							Eigen::Vector3d v1 = edge_list[index_k].v.row(0), v2 = edge_list[index_k].v.row(1);
							if (isEquation(joint_list[i].center, v1))
								ek = v2;
							else if (isEquation(joint_list[i].center, v2))
								ek = v1;
							else
								std::cout << "Error in getTopology() calculate is_same_side" << std::endl;

							double dis = pointToTriangleDistance(p, joint_list[i].center, ej, ek);
							dis = dis - std::max(edge_list[index_j].r_outer, edge_list[index_k].r_outer);
							if (dis < minDis)
								minDis = dis;

							/*std::pair<double, double> dis_pair = pointToPartialPlaneDistance(p, joint_list[i].center, ej, ek);
							double dis1 = dis_pair.first - std::max(edge_list[index_j].r_outer, edge_list[index_k].r_outer);
							double dis2 = dis_pair.second - std::max(edge_list[index_j].r_outer, edge_list[index_k].r_outer);
							if (dis1 < minDis1)
								minDis1 = dis1;
							if (dis2 < minDis2)
								minDis2 = dis2;*/
						}
					}

					if (minDis >= 0)
					{
						result_joint = std::max(result_joint, minDis);
					}
					else if (result_joint < 0)
						result_joint = result_joint * (-1) * tanh(minDis);

				}
			}
			result.push_back(result_joint);
		}

		double sum_KS = 0;	//总的KS函数值
		double sum_smf = 0.1;
		for (int i = 0; i < result.size(); i++)
		{
			sum_KS += exp(-1 / sum_smf * result[i]);
		}
		sum_KS = -sum_smf * log(sum_KS);
		return sum_KS;
	};


	Eigen::VectorXd B(GV.rows());
#pragma omp parallel for
	for (int i = 0; i < GV.rows(); i++) {
		double indicator = function(GV(i, 0), GV(i, 1), GV(i, 2));
		B(i) = indicator;
	}
	igl::copyleft::marching_cubes(B, GV, (res_mesh(0) + 2), (res_mesh(1) + 2), (res_mesh(2) + 2), V, F);

	//平移min到原点
	if (V.rows() == 0 || !V.allFinite())
		throw std::runtime_error("No finite isosurface; check the full-cell skeleton, radius and mesh resolution");
	Eigen::Vector3d minV = V.colwise().minCoeff();
	Eigen::Vector3d maxV = V.colwise().maxCoeff();

	//Eigen::Vector3d move = minV + (maxV - minV) / 2.0;
	Eigen::Vector3d move = minV;
	for (int j = 0; j < V.rows(); j++)
	{
		V(j, 0) -= move(0);
		V(j, 1) -= move(1);
		V(j, 2) -= move(2);
	}

	//缩放为[0,0,0]开始的边长为1的立方体
	Eigen::Vector3d factor = maxV - minV;
	//double xfac = ((double)size_x / factor(0));
	//double yfac = ((double)size_y / factor(1));
	//double zfac = ((double)size_z / factor(2));
	double xfac = 50.0 / factor(0), yfac = 50.0 / factor(1), zfac = 50.0 / factor(2);
	for (int j = 0; j < V.rows(); j++)
	{

		V(j, 0) *= xfac;
		V(j, 1) *= yfac;
		V(j, 2) *= zfac;
	}
}
