#include"Parameterization.h"
int main()
{
	MyMesh inputMesh;
	if (!OpenMesh::IO::read_mesh(inputMesh, "D:/data/Parameterization/nefertiti.off")) {
		std::cerr << "Error: Cannot read mesh from file." << std::endl;
		return 1;
	}
	Parameterization a(inputMesh);
	a.apply();
}
