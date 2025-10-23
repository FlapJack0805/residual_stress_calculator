#include <CGAL/Simple_cartesian.h>
#include <memory>
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

#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
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

#include "nlopt.hpp"

#include <filesystem>
#include <string>

using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
//using Kernel = CGAL::Simple_cartesian<double>;
using Domain = CGAL::Polyhedral_mesh_domain_with_features_3<Kernel, CGAL::Surface_mesh<Kernel::Point_3>>;
using Tr = CGAL::Mesh_triangulation_3<Domain>::type;
using C3t3 = CGAL::Mesh_complex_3_in_triangulation_3<Tr>;
using Criteria = CGAL::Mesh_criteria_3<Tr>;
using Point = Kernel::Point_3;
using Vector = Kernel::Vector_3;
using Surface_mesh = CGAL::Surface_mesh<Point>;
using Vertex_descriptor = boost::graph_traits<Surface_mesh>::vertex_descriptor;
using Primitive = CGAL::AABB_face_graph_triangle_primitive<Surface_mesh>;
using Traits = CGAL::AABB_traits_3<Kernel, Primitive>;
using Tree = CGAL::AABB_tree<Traits>;

#pragma once
class FEA; 
class MeshHandler;


/*
 * UN_CUT means the node hasn't been cut yet
 * CUT means the node has already been cut out
 * CUT_NEXT means the node is about to be cut out
*/
enum class NodeStatus 
{
	UN_CUT = 0,
	CUT_NEXT,
	CUT,
	FORCE_NODE
};


struct MyVector
{
	double x, y, z;
};


/*
 * only stores normal stresses
 * assume sheer stresses have a negligible impact
*/
struct StressTensor
{
	double xx, yy, zz;
};


struct Node
{
	Tr::Vertex_handle handle;
	Point position;
	MyVector normal_vector;
	StressTensor stress;
	NodeStatus status;
	bool is_boundary;
	std::vector<MyVector> deformations; // store all deformations so we can backtrack at the end to get to origional location

	Node(Tr::Vertex_handle handle_in, Point position_in, MyVector normal_vector_in)
	: handle(handle_in), position(position_in), normal_vector(normal_vector_in), status(NodeStatus::UN_CUT), is_boundary(false) {}

	Node() : status(NodeStatus::UN_CUT), is_boundary(false) {}
};


struct OptimizationContext 
{
    MeshHandler* handler;
    Surface_mesh* goal_mesh;
};



class MeshHandler
{
	C3t3 volume_mesh;
	std::vector<std::shared_ptr<Node>> nodes; // each index is the node id
	void set_tool_path(const std::filesystem::path& tool_path_stl);

public:
	MeshHandler(const std::filesystem::path& mesh_stl, const std::filesystem::path& boundary_stl);
	void set_stress_in_cut(const std::filesystem::path& tool_path_stl, const std::filesystem::path& deformed_mesh);
	void print_summary() const;

	inline const std::vector<std::shared_ptr<Node>>& get_nodes() const { return nodes; }
	inline std::shared_ptr<Node> get_node(size_t id) { return nodes.at(id); }
	inline size_t get_node_count() const { return nodes.size(); }
	inline const C3t3& get_c3t3() const { return volume_mesh; }
	inline C3t3& get_c3t3() { return volume_mesh; }
};
