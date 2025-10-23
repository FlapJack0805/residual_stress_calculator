#include "mesh_handler.hpp"
#include "fea_solver.hpp"




/*
 * This is the constructor for the MeshHandler class.
 * It takes in an stl file for the material mesh and for the boundary mesh which then makes our volumetric mesh that we will use all of our computations for 
 * as well as marking every node that is a boundary node.
*/
MeshHandler::MeshHandler(const std::filesystem::path& mesh_stl, const std::filesystem::path& boundary_stl)
{
    Surface_mesh surface_mesh;
    Surface_mesh boundary_mesh;

    // --- Load the main and boundary meshes ---
    if (!CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(mesh_stl.string(), surface_mesh)) {
        throw std::runtime_error("Failed to load mesh: " + mesh_stl.string());
    }

    if (!CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(boundary_stl.string(), boundary_mesh)) {
        throw std::runtime_error("Failed to load boundary: " + boundary_stl.string());
    }

    if (!CGAL::Polygon_mesh_processing::is_outward_oriented(boundary_mesh)) {
        CGAL::Polygon_mesh_processing::reverse_face_orientations(boundary_mesh);
    }

    // --- Build the volume mesh (tetrahedralization) ---
    Domain domain(surface_mesh);
    Criteria criteria(CGAL::parameters::edge_size(5)
                      .facet_angle(5)
                      .facet_distance(1)
                      .cell_radius_edge_ratio(3.0)
                      .cell_size(10));

    volume_mesh = CGAL::make_mesh_3<C3t3>(domain, criteria);
    std::cout << "[MeshHandler] Volume mesh created: "
              << volume_mesh.triangulation().number_of_vertices()
              << " vertices\n";

    // --- Build point-in-mesh tester for the boundary ---
    using Primitive = CGAL::AABB_face_graph_triangle_primitive<Surface_mesh>;
    using Traits    = CGAL::AABB_traits_3<Kernel, Primitive>;
    using Tree      = CGAL::AABB_tree<Traits>;

    Tree boundary_tree(faces(boundary_mesh).begin(), faces(boundary_mesh).end(), boundary_mesh);
    boundary_tree.accelerate_distance_queries();

    CGAL::Side_of_triangle_mesh<Surface_mesh, Kernel> inside_test(boundary_tree);

    // --- Iterate through all finite vertices in the tetrahedral mesh ---
    const auto& triangulation = volume_mesh.triangulation();
    size_t node_id = 0;
    size_t num_boundary_nodes = 0; // This is purely for testing
    for (auto vit = triangulation.finite_vertices_begin();
         vit != triangulation.finite_vertices_end(); ++vit, ++node_id)
    {
        Point p = vit->point().point();
        MyVector n = {0, 0, 1}; // placeholder (you can compute vertex normal later)
        auto node = std::make_shared<Node>(vit, p, n);
        nodes.push_back(node);

        // --- Classify boundary membership ---
        if (inside_test(p) == CGAL::ON_BOUNDARY ||
            inside_test(p) == CGAL::ON_BOUNDED_SIDE)
        {
            node->is_boundary = true;
            num_boundary_nodes++;
        }
    }

    std::cout << "[MeshHandler] Total nodes: " << nodes.size() << "\n";
    std::cout << "[MeshHandler] Boundary nodes: " << num_boundary_nodes << "\n";
}



/*
 * Takes in a tool path and updates the status of all of the nodes that are either going to be cut or used as force nodes.
 * Everything inside of the tool path mesh will be cut and everything that is on the boundary of it will be used as force nodes.
*/
void MeshHandler::set_tool_path(const std::filesystem::path& tool_path_stl)
{
    Surface_mesh tool_path_mesh;
    if (!CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(tool_path_stl.string(), tool_path_mesh)) 
    {
        throw std::runtime_error("Failed to load mesh: " + tool_path_stl.string());
    }


    Tree tree(faces(tool_path_mesh).first, faces(tool_path_mesh).second, tool_path_mesh);
    tree.accelerate_distance_queries();

    CGAL::Side_of_triangle_mesh<Surface_mesh, Kernel> inside_test(tree);


    for (auto& node : nodes) {
        Point p = node->position;

        CGAL::Bounded_side side = inside_test(p);

        if (side == CGAL::ON_BOUNDED_SIDE) {
            node->status = NodeStatus::CUT_NEXT;
        } 
        else if (side == CGAL::ON_BOUNDARY) {
            node->status = NodeStatus::FORCE_NODE;
        }
    }
}


static double intersection_volume(const C3t3& c3,
                           const Surface_mesh& shell)
{
    // 1. Build an AABB tree and a point-inside tester on the shell
    using Primitive = CGAL::AABB_face_graph_triangle_primitive<Surface_mesh>;
    using Traits    = CGAL::AABB_traits_3<Kernel, Primitive>;
    using Tree      = CGAL::AABB_tree<Traits>;

	CGAL::AABB_tree<CGAL::AABB_traits_3<Kernel, CGAL::AABB_face_graph_triangle_primitive<Surface_mesh>>> 
	tree(faces(shell).begin(), faces(shell).end(), shell);
	// Construct side tester
	CGAL::Side_of_triangle_mesh<Surface_mesh, Kernel> inside_test(tree);

    // 2. Sweep all cells in the complex and accumulate volume
    Kernel::FT vol = 0;
    for(const auto& cell : c3.cells_in_complex())
    {
        // barycenter of the tet
        const Kernel::Point_3& p0 = cell->vertex(0)->point().point();
        const Kernel::Point_3& p1 = cell->vertex(1)->point().point();
        const Kernel::Point_3& p2 = cell->vertex(2)->point().point();
        const Kernel::Point_3& p3 = cell->vertex(3)->point().point();
        Kernel::Point_3 bary(
            CGAL::centroid(p0, p1, p2, p3));

        if( inside_test(bary) == CGAL::ON_BOUNDED_SIDE )
            vol += CGAL::volume(p0, p1, p2, p3);
    }
    return CGAL::to_double(vol);     // intersection volume
}


static double objective_func(const std::vector<double>& x, std::vector<double>& grad, void* data)
{
    OptimizationContext* context = static_cast<OptimizationContext*>(data);
    MeshHandler* handler = context->handler;
    Surface_mesh* goal_mesh = context->goal_mesh;

    double pressure = x[0];

    FEA fea(*handler);
    C3t3 deformed_mesh = fea.apply_trench_force(pressure);

    double volume = intersection_volume(deformed_mesh, *goal_mesh);
    std::cout << "[NLopt] pressure = " << pressure << " MPa → volume = " << volume << "\n";
    return volume; // NLopt will maximize this
}


void MeshHandler::set_stress_in_cut(const std::filesystem::path& tool_path_stl, const std::filesystem::path& deformed_stl)
{
    set_tool_path(tool_path_stl); // sets all values in nodes vector to what we need to simulate the cut

    Surface_mesh deformed_mesh;

    if (!CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(deformed_stl.string(), deformed_mesh)) 
    {
        throw std::runtime_error("Failed to load mesh: " + deformed_stl.string());
    }

    OptimizationContext context{this, &deformed_mesh};

    nlopt::opt opt(nlopt::LN_COBYLA, 1);  // 1 parameter: pressure
    opt.set_max_objective(objective_func, &context);
    opt.set_lower_bounds(-5000.0);
    opt.set_upper_bounds(5000.0);
    opt.set_xtol_rel(1e-3);

    std::vector<double> x = {100.0};
    double max_volume;
    opt.optimize(x, max_volume);

    std::cout << "[Optimization] Best pressure = "
              << x[0] << " MPa (volume = " << max_volume << ")\n";


}


/*
 * Prints a summary of the entire volumetric mesh
*/
void MeshHandler::print_summary() const
{
    size_t total_nodes = nodes.size();
    size_t boundary_count = 0;
    size_t force_count = 0;
    size_t uncut_count = 0;
    size_t cutnext_count = 0;
    size_t cut_count = 0;

    for (const auto& node : nodes) {
        if (node->is_boundary) boundary_count++;

        switch (node->status) {
            case NodeStatus::UN_CUT:    uncut_count++; break;
            case NodeStatus::CUT_NEXT:  cutnext_count++; break;
            case NodeStatus::CUT:       cut_count++; break;
            case NodeStatus::FORCE_NODE: force_count++; break;
        }
    }

    auto pct = [&](size_t n) -> double {
        return total_nodes ? (100.0 * n / total_nodes) : 0.0;
    };

    std::cout << "\n========== Mesh Summary ==========\n";
    std::cout << "Total nodes:      " << total_nodes << "\n";
    std::cout << "Boundary nodes:   " << boundary_count << " (" << pct(boundary_count) << "%)\n";
    std::cout << "Force nodes:      " << force_count << " (" << pct(force_count) << "%)\n";
    std::cout << "Status breakdown:\n";
    std::cout << "   UN_CUT:        " << uncut_count   << " (" << pct(uncut_count)   << "%)\n";
    std::cout << "   CUT_NEXT:      " << cutnext_count << " (" << pct(cutnext_count) << "%)\n";
    std::cout << "   CUT:           " << cut_count     << " (" << pct(cut_count)     << "%)\n";
    std::cout << "==================================\n\n";
}
