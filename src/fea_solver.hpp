#include "mesh_handler.hpp"

class FEA;

class FEA
{
	//later we will likely have different materials we can work with so we will customize it here but for simplicity now I will ignore it
	MeshHandler& mesh;
	
public:
	FEA(MeshHandler& mesh_in);
	C3t3 apply_trench_force(double pressure);
	void get_cut_stress_field(double pressure);
	void get_full_stress_field();
};
