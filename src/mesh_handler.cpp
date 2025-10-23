#include "mesh_handler.hpp"




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
        UnitVector n = {0, 0, 1}; // placeholder (you can compute vertex normal later)
        auto node = std::make_shared<Node>(Vertex_descriptor(), p, n);
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
