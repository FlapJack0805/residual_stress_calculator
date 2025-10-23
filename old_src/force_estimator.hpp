#include "CLI/CLI.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <regex>
#include <cmath>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Search_traits_2.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_iso_box.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/Surface_mesh.h>
#include "CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h"
#include <CGAL/Polygon_mesh_processing/transform.h>
#include "CGAL/Polygon_mesh_processing/orientation.h"
#include "CGAL/Polygon_mesh_processing/self_intersections.h"
#include <CGAL/K_neighbor_search.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>

#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>   // <- use this
#include <CGAL/Mesh_triangulation_3.h>               
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>  
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Mesh_criteria_3.h>                    
#include <CGAL/make_mesh_3.h>                        
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/aff_transformation_tags.h>

#include "CGAL/boost/graph/IO/STL.h"
#include "CGAL/IO/io.h"

#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/boost/graph/helpers.h> // for vertices(), etc.
#include <vector>
#include <unordered_map>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Polygon_mesh_processing/transform.h>
#include <CGAL/Kernel/hash_functions.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Kernel/global_functions_3.h>
#include <CGAL/perturb_mesh_3.h>   

#include <boost/property_map/property_map.hpp>

#include <Eigen/Dense>

#include "nlopt.hpp"


using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
//using Kernel = CGAL::Simple_cartesian<double>;
using Domain = CGAL::Polyhedral_mesh_domain_with_features_3<Kernel, CGAL::Surface_mesh<Kernel::Point_3>>;
using Tr = CGAL::Mesh_triangulation_3<Domain>::type;
using C3t3 = CGAL::Mesh_complex_3_in_triangulation_3<Tr>;
using Criteria = CGAL::Mesh_criteria_3<Tr>;
using Point = Kernel::Point_3;
using Vector = Kernel::Vector_3;
using Traits = CGAL::Search_traits_2<Kernel>; using KD_tree =  CGAL::Kd_tree<Traits>;
using Search = CGAL::Orthogonal_k_neighbor_search<Traits>;
using Surface_mesh = CGAL::Surface_mesh<Point>;
using Vertex_descriptor = boost::graph_traits<Surface_mesh>::vertex_descriptor;


struct UnitVector
{
	double x, y;
};


struct StressTensor
{
	double xx, yy, zz;
};


/*
 * UN_CUT means the node hasn't been cut yet
 * CUT means the node has already been cut out
 * CUT_NEXT means the node is about to be cut out
*/
enum class NodeStatus 
{
	UN_CUT = 0,
	CUT_NEXT,
	CUT
};


struct Node
{
	Vertex_descriptor v;
	Point position;
	UnitVector normal_vector;
	StressTensor stress;
	NodeStatus is_cut;

	Node(Vertex_descriptor v_in, Point position_in, UnitVector normal_vector_in)
	: v(v_in), position(position_in), normal_vector(normal_vector_in), is_cut(NodeStatus::UN_CUT) {}

	Node() : is_cut(NodeStatus::UN_CUT) {}
};


struct Edge 
{
	size_t a, b;
	bool operator==(const Edge o) const
	{
		return a == o.a && b == o.b;
	}
};


struct EdgeHash
{
	size_t operator()(const Edge & e) const 
	{
		return (e.a*1315423911u) ^ e.b;
	}
};


struct Ten 
{ 
	size_t v[10]; 
};

class ToolPath
{
	KD_tree tool_path_tree;
	Surface_mesh tool_path_mesh;
	Point start_point;

public:
	ToolPath(const std::filesystem::path& stl_file);
	void create_nodes(C3t3& mesh, std::unordered_map<Tr::Vertex_handle, size_t> &node_id_map);
};


class ResidualCalculator 
{
	std::vector<Node> force_nodes;
	std::vector<std::pair<std::size_t,int>> trench_faces;
	std::vector<std::array<std::size_t,4>> elem_connect;
	std::vector<Tr::Vertex_handle> boundary_nodes;
	std::unordered_map<Edge, size_t, EdgeHash> edge_mid_map;
	std::unordered_map<Tr::Vertex_handle, size_t> node_id_map;
	std::unordered_map<size_t, Tr::Vertex_handle> id_node_map;
	std::unordered_map<size_t, Point> id_point_map;
	std::unordered_map<Point, Tr::Vertex_handle> point_to_vertex_handle;
	std::vector<Ten> elem10;
	C3t3 mesh_c3t3;
	Surface_mesh goal_mesh;
	std::unordered_map<Tr::Vertex_handle, StressTensor> get_stress_tensor();
	std::vector<size_t> get_elems_in_mesh(const Surface_mesh &tool_path);


public:
	ResidualCalculator(const std::filesystem::path& object_stl, const std::filesystem::path& tool_path_stl, 
		    const std::filesystem::path& boundary_stl, const std::filesystem::path& goal_stl);
	std::unordered_map<Tr::Vertex_handle, StressTensor> stress_estimator();
	double obj(double force, double t_x, double t_y, double t_z,
				       double r_x, double r_y, double r_z);
	Surface_mesh simulate_cut(const Surface_mesh &tool_path);
	std::unordered_map<Tr::Vertex_handle, StressTensor> full_material_stress_estimator();
};
