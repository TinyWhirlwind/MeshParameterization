#ifndef PARAMETERIZATION_H
#define PARAMETERIZATION_H
#include <iostream>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Mesh/Traits.hh>
#include <Eigen/Dense>
typedef OpenMesh::TriMesh_ArrayKernelT<>  MyMesh;
class Parameterization
{
	enum ParamType
	{
		Param_BM,
		Param_LSCM
	};
public:
	Parameterization();
	Parameterization(MyMesh mesh);
	~Parameterization();
public:
	void apply(ParamType type);
	void TextureMapping();
	MyMesh getMesh();
private:
	class PImpl;
	std::shared_ptr<PImpl> impl_;
};
#endif