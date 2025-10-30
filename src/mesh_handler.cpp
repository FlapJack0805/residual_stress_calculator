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
    Criteria criteria(CGAL::parameters::edge_size(3)
                      .facet_angle(20)
                      .facet_distance(0.25)
                      .cell_radius_edge_ratio(2.0)
                      .cell_size(5));

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


    for (auto& node : nodes) 
    {
        Point p = node->position;

        CGAL::Bounded_side side = inside_test(p);

        if (side == CGAL::ON_BOUNDED_SIDE) 
        {
            node->status = NodeStatus::CUT_NEXT;
        } 
        else if (side == CGAL::ON_BOUNDARY) 
        {
            if (side == CGAL::ON_BOUNDARY)
            {
                node->status = NodeStatus::FORCE_NODE;
            }
        }
    }
}



/*
 * Custom function to get the intersection volume between a volumetric mesh and a surface mesh because CGAL doesn't provide one
 * Works by testing for every element in the volumetric mesh to see if it's in the surface mesh and sum up all the elements that are
*/
static std::pair<double, double> intersection_volume(const C3t3& c3,
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
    Kernel::FT intersection_volume = 0;
    Kernel::FT total_volume = 0;
    for(const auto& cell : c3.cells_in_complex())
    {
        // barycenter of the tet
        const Kernel::Point_3& p0 = cell->vertex(0)->point().point();
        const Kernel::Point_3& p1 = cell->vertex(1)->point().point();
        const Kernel::Point_3& p2 = cell->vertex(2)->point().point();
        const Kernel::Point_3& p3 = cell->vertex(3)->point().point();
        Kernel::Point_3 bary(
            CGAL::centroid(p0, p1, p2, p3));
        Kernel::FT volume = CGAL::volume(p0, p1, p2, p3);
        total_volume += volume;
        if( inside_test(bary) == CGAL::ON_BOUNDED_SIDE )
        {
            intersection_volume += volume;
        }
    }
    return {CGAL::to_double(intersection_volume), CGAL::to_double(total_volume)};     // intersection volume
    
    /*
    namespace PMP = CGAL::Polygon_mesh_processing;

    double vol = 0.0;
    auto& tr = c3.triangulation();

    for (auto cit = tr.finite_cells_begin(); cit != tr.finite_cells_end(); ++cit)
    {
        Surface_mesh tet_mesh;
        CGAL::make_tetrahedron(
            cit->vertex(0)->point().point(),
            cit->vertex(1)->point().point(),
            cit->vertex(2)->point().point(),
            cit->vertex(3)->point().point(),
            tet_mesh);

        // Copy since clipping modifies input
        Surface_mesh clipped = tet_mesh;
        Surface_mesh shell_copy = shell;

        PMP::clip(clipped, shell_copy, PMP::parameters::clip_volume(true));

        if (!clipped.is_empty())
            vol += PMP::volume(clipped);
    }

    return vol;
    */
}


/*
 * Objection function used by the optimizer
 * It applies a given pressure to the mesh and outputs the intersection volume of the new deformed mesh with the goal mesh
 * with a penalty applied if for not matching the goal mesh's volume
*/
static double objective_func(const std::vector<double>& x, std::vector<double>& grad, void* data)
{
    OptimizationContext* context = static_cast<OptimizationContext*>(data);
    MeshHandler* handler = context->handler;
    Surface_mesh* goal_mesh = context->goal_mesh;
    double goal_volume = context->goal_volume;

    double pressure = x[0];

    FEA fea(*handler);
    C3t3 deformed_mesh = fea.apply_trench_force(pressure);

    std::pair<double, double> volumes = intersection_volume(deformed_mesh, *goal_mesh);
    double intersection_volume = volumes.first;
    double deformed_volume = volumes.second;
    double volume_mismatch = std::abs(goal_volume - deformed_volume);
    double objective = intersection_volume - context->lambda * volume_mismatch;
    std::cout << "[OBJ] pressure=" << pressure
              << " | overlap=" << intersection_volume
              << " | dVol=" << deformed_volume
              << " | penalty=" << context->lambda * volume_mismatch
              << " | total=" << objective 
              << "\n";
    return objective; // NLopt will maximize this

}


/*
 * This function takes in a tool path mesh of the cut and the deformed mesh after the cut and calculates the stress
 * It stores the stress values in this class so it doesn't return anything
*/
void MeshHandler::set_stress_in_cut(const std::filesystem::path& tool_path_stl, const std::filesystem::path& deformed_stl)
{
    export_surface_mesh("starting_mesh.stl");
    set_tool_path(tool_path_stl); // sets all values in nodes vector to what we need to simulate the cut
    
    std::cout << "Summary after tool path" << std::endl;
    print_summary();

    Surface_mesh deformed_mesh;

    if (!CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(deformed_stl.string(), deformed_mesh)) 
    {
        throw std::runtime_error("Failed to load mesh: " + deformed_stl.string());
    }

    OptimizationContext context{this, &deformed_mesh, CGAL::Polygon_mesh_processing::volume(deformed_mesh), 1};

    nlopt::opt opt(nlopt::LN_COBYLA, 1);  // 1 parameter: pressure
    opt.set_max_objective(objective_func, &context);
    opt.set_lower_bounds(-5000.0);
    opt.set_upper_bounds(5000.0);
    opt.set_xtol_rel(1e-3);

    std::vector<double> x = {100.0};
    double max_volume;
    opt.optimize(x, max_volume);

    for (std::shared_ptr<Node> node : nodes)
    {
        if (node->status == NodeStatus::CUT_NEXT)
        {
            node->status = NodeStatus::CUT;
        }

    }

    std::cout << "[Optimization] Best pressure = "
              << x[0] << " MPa (volume = " << max_volume << ")\n";


    export_surface_mesh("debug_mesh.stl");

    FEA fea(*this);

    //This doesn't work yet
    //fea.get_cut_stress_field(x[0]);

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

    for (const auto& node : nodes) 
    {
        if (node->is_boundary) boundary_count++;

        switch (node->status) 
        {
            case NodeStatus::UN_CUT:    uncut_count++; break;
            case NodeStatus::CUT_NEXT:  cutnext_count++; break;
            case NodeStatus::CUT:       cut_count++; break;
            case NodeStatus::FORCE_NODE: force_count++; break;
        }
    }

    auto pct = [&](size_t n) -> double 
    {
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


/*
 * Outputs the current volumetric mesh stored in the class as an stl file
 * This is very helpful for testing
*/
void MeshHandler::export_surface_mesh(const std::filesystem::path& out_path) const
{
    Surface_mesh surface_mesh;
    CGAL::facets_in_complex_3_to_triangle_mesh(volume_mesh, surface_mesh);

    if (!CGAL::IO::write_polygon_mesh(out_path.string(), surface_mesh, CGAL::parameters::stream_precision(17))) {
        throw std::runtime_error("Failed to write surface mesh to " + out_path.string());
    }

    std::cout << "[MeshHandler] Exported surface mesh to " << out_path << std::endl;
}
