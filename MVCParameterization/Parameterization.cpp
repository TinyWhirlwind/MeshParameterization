#include "Parameterization.h"
#include "Eigen/Core"
#include "Eigen/Sparse"
typedef Eigen::SparseMatrix<double> SparseMat;
typedef Eigen::Triplet<double> T;
class Parameterization::PImpl
{
public:
	PImpl()
	{
	}
	PImpl(MyMesh mesh)
	{
		mesh_ = mesh;
	}
	~PImpl() {}

public:
	void TutteEmbedding(std::map<int, MyMesh::Point>& map_list);
	void CheckBoundarys(std::vector<int>& boundaryPts,float& boundary_length);
	void MapBoundaryToCircle(std::vector<int> boundaryPts, float R, std::map<int,MyMesh::Point>& result);
	void LinkMesh(std::map<int, MyMesh::Point> map_list);
public:
	MyMesh mesh_;
	MyMesh mesh_2d_;
};
Parameterization::Parameterization()
{
	impl_.reset(new PImpl());
}

Parameterization::Parameterization(MyMesh mesh)
{
	impl_.reset(new PImpl(mesh));
}

Parameterization::~Parameterization(){}


void Parameterization::apply(ParamType type)
{
	std::vector<int> boundaryPts;
	float boundary_length = 0.0;
	impl_->CheckBoundarys(boundaryPts, boundary_length);
	std::map<int, MyMesh::Point> map_result;
	switch (type)
	{
	case Param_BM:
		float R = boundary_length / (2 * M_PI);
		impl_->MapBoundaryToCircle(boundaryPts, R, map_result);
		impl_->TutteEmbedding(map_result);
	case Param_LSCM:

	}
	
	impl_->LinkMesh(map_result);
}

void Parameterization::TextureMapping()
{
	OpenMesh::VPropHandleT<MyMesh::TexCoord3D> texcoords;
	OpenMesh::FPropHandleT<MyMesh::Normal> normals;
	if (!impl_->mesh_.get_property_handle(texcoords, "texcoords"))
	{
		impl_->mesh_.add_property(texcoords, "texcoords");
	}
	if (!impl_->mesh_.get_property_handle(normals, "normals"))
	{
		impl_->mesh_.add_property(normals, "normals");
	}
	//设置纹理坐标
	for (auto v_it : impl_->mesh_.vertices())
	{
		impl_->mesh_.property(texcoords, v_it) = MyMesh::TexCoord3D(0.0, 0.0, 0.0);
	}
	for (auto f_it : impl_->mesh_.faces())
	{
		impl_->mesh_.property(normals, f_it) = MyMesh::Normal(0.0, 0.0, 0.0);
	}
}

MyMesh Parameterization::getMesh()
{
	return impl_->mesh_2d_;
}

void Parameterization::PImpl::CheckBoundarys(std::vector<int>& boundaryPts, float& boundary_length)
{
	std::vector<MyMesh::HalfedgeHandle> half_list;
	for (MyMesh::HalfedgeIter h_it = mesh_.halfedges_begin(); h_it != mesh_.halfedges_end(); ++h_it)
	{
		if (h_it->is_boundary())
		{
			MyMesh::HalfedgeHandle h = *h_it;
			do {
				half_list.push_back(h);
				const MyMesh::Point& vTo = mesh_.point(mesh_.to_vertex_handle(h));
				const MyMesh::Point& vFrom = mesh_.point(mesh_.from_vertex_handle(h));
				h = mesh_.next_halfedge_handle(h);
			} while (h != *h_it);
			break;
		}
	}
	boundary_length = 0.0f;
	for (auto h_it : half_list)
	{
		MyMesh::VertexHandle vFrom = mesh_.from_vertex_handle(h_it);
		boundaryPts.push_back(vFrom.idx());
		boundary_length += mesh_.calc_edge_length(h_it);
	}
}

void Parameterization::PImpl::MapBoundaryToCircle(std::vector<int> boundaryPts, float R, std::map<int, MyMesh::Point>& result)
{
	result.clear();
	int boundary_num = boundaryPts.size();
	float part = 2 * M_PI / boundary_num;
	for (int i = 0; i < boundary_num; i++)
	{
		MyMesh::Point circleP;
		circleP[0] = R * cos(i * part);
		circleP[1] = R * sin(i * part);
		circleP[2] = 0.0f;
		result[boundaryPts[i]] = circleP;
	}
}

void Parameterization::PImpl::TutteEmbedding(std::map<int, MyMesh::Point>& map_list)
{
	int n = mesh_.n_vertices();
	SparseMat coefficientMat(n,n);//系数矩阵
	std::vector<T> triplets;
	
	Eigen::MatrixXd targetVector(n, 2);
	targetVector.setZero();
	Eigen::MatrixXd resultVector(n, 2);
	for (MyMesh::VertexIter v_it = mesh_.vertices_begin(); v_it != mesh_.vertices_end(); ++v_it)
	{
		int  i = v_it->idx();
		if (mesh_.is_boundary(*v_it))
		{
			targetVector(i, 0) = map_list[i][0];
			targetVector(i, 1) = map_list[i][1];
			triplets.push_back(T(i, i, 1));
		}
		else
		{
			int v_degree = 0;
			for (MyMesh::VertexVertexIter vv_it = mesh_.vv_begin(*v_it); vv_it .is_valid(); ++vv_it) 
			{
				int j = vv_it->idx();
				++v_degree;
				triplets.push_back(T(i, j, 1.0));
			}
			triplets.push_back(T(i, i, -double(v_degree)));
		}
	}
	coefficientMat.resize(n, n);
	coefficientMat.setFromTriplets(triplets.begin(), triplets.end());
	Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> solver;
	solver.compute(coefficientMat);
	resultVector = solver.solve(targetVector);
	for (int i = 0; i < n; i++)
	{
		if (!mesh_.is_boundary(mesh_.vertex_handle(i)))
		{
			map_list[i][0] = resultVector.coeff(i, 0);
			map_list[i][1] = resultVector.coeff(i, 1);
			map_list[i][2] = 0.0f;
		}
	}
}
void Parameterization::PImpl::LinkMesh(std::map<int, MyMesh::Point> map_list)
{
	for (auto mp : map_list)
	{
		mesh_2d_.add_vertex(mp.second);
	}
	for (MyMesh::FaceIter f_it = mesh_.faces_begin(); f_it != mesh_.faces_end(); ++f_it)
	{
		auto v_iter = mesh_.fv_begin(*f_it);
		std::vector<MyMesh::VertexHandle> v_list;
		for (int i = 0; i < 3; i++)
		{
			v_list.push_back(*v_iter);
			++v_iter;
		}
		mesh_2d_.add_face(v_list[0], v_list[1], v_list[2]);
	}
}


