#include "force_estimator.hpp"

int main(int argc, char ** argv)
{
	CLI::App app { "tool_path_extractor" };

	//std::string tool_path_nc;
	std::string tool_path_stl;
	std::string clamp_point_stl;
	std::string pre_deformation_stl;
	std::string post_deformation_stl;

	app.add_option("--pre_deformation_stl", pre_deformation_stl, ".stl pre deformation file")->required()->check(CLI::ExistingFile);
	app.add_option("--post_deformation_stl", post_deformation_stl, ".stl post deformation file")->required()->check(CLI::ExistingFile);
	app.add_option("--tool_path_stl", tool_path_stl, ".stl tool path file")->required()->check(CLI::ExistingFile);
	app.add_option("--clamp_points_mesh_stl", clamp_point_stl, ".stl clamp point file")->required()->check(CLI::ExistingFile);
                                                     //"Format: \"(x,y,z)\"")->required();

	CLI11_PARSE(app, argc, argv);

	ResidualCalculator residual_calculator(pre_deformation_stl, tool_path_stl, clamp_point_stl, post_deformation_stl);
	std::unordered_map<Tr::Vertex_handle, StressTensor> trench_stresses = residual_calculator.stress_estimator();

	std::cout << "stress map size: " << trench_stresses.size() << '\n';
	for (const std::pair<Tr::Vertex_handle, StressTensor> &trench_stress: trench_stresses)
	{
		std::cout << "nxx: " << trench_stress.second.xx
			<< " nyy: " << trench_stress.second.yy
			<< " nzz: " << trench_stress.second.zz
			<< '\n';
	}
}


