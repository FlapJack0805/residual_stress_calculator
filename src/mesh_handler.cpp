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

    if (!CGAL::is_closed(surface_mesh)) 
    {
        std::cerr << "SURFACE MESH IS NOT CLOSED!" << std::endl;
    }

    // --- Build the volume mesh (tetrahedralization) ---
    Domain domain(surface_mesh); 
    domain.detect_features();

    Criteria criteria(
        CGAL::parameters::edge_size(2)
        .facet_angle(20)
        .facet_distance(0)     // do not move boundary
        .facet_size(0)         // do not change triangle density
        .cell_radius_edge_ratio(2.0)
        .cell_size(4)
    );

    /*
    std::vector<std::tuple<Kernel::Weighted_point_3, int, Domain::Surface_patch_index>> surface_mesh_points;
    for (auto v : surface_mesh.vertices())
    {
         Kernel::Weighted_point_3 wp(surface_mesh.point(v));
        surface_mesh_points.emplace_back(wp, 2 , Domain::Surface_patch_index(0));
    }
    */

    volume_mesh = CGAL::make_mesh_3<C3t3>(
        domain, 
        criteria
    );
    
    std::cout << "[MeshHandler] Volume mesh created: "
              << volume_mesh.triangulation().number_of_vertices()
              << " vertices\n";

    export_surface_mesh("mesh_before_optimization.stl", &volume_mesh);

    // --- Build point-in-mesh tester for the boundary ---
    using Primitive = CGAL::AABB_face_graph_triangle_primitive<Surface_mesh>;
    using Traits    = CGAL::AABB_traits_3<Kernel, Primitive>;
    using Tree      = CGAL::AABB_tree<Traits>;

    /*
    Tree sm_tree(faces(surface_mesh).begin(), faces(surface_mesh).end(), surface_mesh);
    sm_tree.accelerate_distance_queries();

    std::unordered_set<Tr::Vertex_handle> boundary_vertices;

    for (auto fit = volume_mesh.facets_in_complex_begin();
         fit != volume_mesh.facets_in_complex_end(); ++fit)
    {
        Tr::Cell_handle c = fit->first;
        int i = fit->second;

        for (int v = 0; v < 4; ++v)
        {
            if (v == i) continue;
            boundary_vertices.insert(c->vertex(v));
        }
    }

    for (auto vh : boundary_vertices)
    {
        auto wp = vh->point();         // weighted point
        Point p = wp.point();          // unweighted
        double w = wp.weight();        // preserve weight

        Point closest = sm_tree.closest_point(p);

        vh->set_point(Kernel::Weighted_point_3(closest, w));
    }
    */

    // NOTE: The boundary nodes refer to the nodes that are right on the clamp and don't move during FEA simulation and exterior nodes are on the exterior

    Tree boundary_tree(faces(boundary_mesh).begin(), faces(boundary_mesh).end(), boundary_mesh);
    boundary_tree.accelerate_distance_queries();

    CGAL::Side_of_triangle_mesh<Surface_mesh, Kernel> inside_test(boundary_tree);

    Tree exterior_tree(faces(surface_mesh).begin(), faces(surface_mesh).end(), surface_mesh);
    exterior_tree.accelerate_distance_queries();

    CGAL::Side_of_triangle_mesh<Surface_mesh, Kernel> boundary_test(boundary_tree);
    CGAL::Side_of_triangle_mesh<Surface_mesh, Kernel> exterior_test(exterior_tree);

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
        if (boundary_test(p) == CGAL::ON_BOUNDARY ||
            boundary_test(p) == CGAL::ON_BOUNDED_SIDE)
        {
            node->is_boundary = true;
            num_boundary_nodes++;
        }

        if (exterior_test(p) == CGAL::ON_BOUNDARY)
        {
            node->is_exterior = true;
        }
    }


    for (auto fit = volume_mesh.facets_in_complex_begin();
         fit != volume_mesh.facets_in_complex_end(); ++fit)
    {
        
    }

    std::cout << "[MeshHandler] Total nodes: " << nodes.size() << "\n";
    std::cout << "[MeshHandler] Boundary nodes: " << num_boundary_nodes << "\n";

}


/*
 * This constructor is for if you already have the meshes made. I really made it so that I could test it with my own boundary nodes and didn't need to build a mesh because 
 * that works weird when working with Abaqus as an absolute struth
*/
MeshHandler::MeshHandler(const Surface_mesh& surface_mesh, const std::vector<Point>& boundary_points)
{
    Domain domain(surface_mesh); 
    domain.detect_features();

    Criteria criteria(
        CGAL::parameters::edge_size(3)
        .facet_angle(20)
        .facet_distance(0)     // do not move boundary
        .facet_size(0)         // do not change triangle density
        .cell_radius_edge_ratio(2.0)
        .cell_size(5)
    );

    /*
    std::vector<std::tuple<Kernel::Weighted_point_3, int, Domain::Surface_patch_index>> surface_mesh_points;
    for (auto v : surface_mesh.vertices())
    {
         Kernel::Weighted_point_3 wp(surface_mesh.point(v));
        surface_mesh_points.emplace_back(wp, 2 , Domain::Surface_patch_index(0));
    }
    */

    volume_mesh = CGAL::make_mesh_3<C3t3>(
        domain, 
        criteria
    );

    std::cout << "[MeshHandler] Volume mesh created: "
              << volume_mesh.triangulation().number_of_vertices()
              << " vertices\n";

    export_surface_mesh("mesh_before_optimization.stl", &volume_mesh);

    // --- Build point-in-mesh tester for the boundary ---
    using Traits        = CGAL::Search_traits_3<Kernel>;
    using Tree          = CGAL::Kd_tree<Traits>;
    using Distance      = CGAL::Orthogonal_k_neighbor_search<Traits>::Distance;
    using Neighbor_search = CGAL::Orthogonal_k_neighbor_search<Traits>;

    /*
    Tree sm_tree(faces(surface_mesh).begin(), faces(surface_mesh).end(), surface_mesh);
    sm_tree.accelerate_distance_queries();

    std::unordered_set<Tr::Vertex_handle> boundary_vertices;

    for (auto fit = volume_mesh.facets_in_complex_begin();
         fit != volume_mesh.facets_in_complex_end(); ++fit)
    {
        Tr::Cell_handle c = fit->first;
        int i = fit->second;

        for (int v = 0; v < 4; ++v)
        {
            if (v == i) continue;
            boundary_vertices.insert(c->vertex(v));
        }
    }

    for (auto vh : boundary_vertices)
    {
        auto wp = vh->point();         // weighted point
        Point p = wp.point();          // unweighted
        double w = wp.weight();        // preserve weight

        Point closest = sm_tree.closest_point(p);

        vh->set_point(Kernel::Weighted_point_3(closest, w));
    }
    */

    // --- Iterate through all finite vertices in the tetrahedral mesh ---
    const auto& triangulation = volume_mesh.triangulation();
    size_t node_id = 0;

    Tree kd_tree(boundary_points.begin(), boundary_points.end());

    size_t num_boundary_nodes = 0;
    double boundary_threshold = 1;

    for (auto vit = triangulation.finite_vertices_begin();
         vit != triangulation.finite_vertices_end();
         ++vit)
    {
        Point p = vit->point().point();

        Neighbor_search search(kd_tree, p, 1);  
        double dist_sq = search.begin()->second;

        auto node = std::make_shared<Node>(vit, p, MyVector{0,0,1});
        nodes.push_back(node);

        if (dist_sq < boundary_threshold * boundary_threshold)
        {
            node->is_boundary = true;
            num_boundary_nodes++;
        }
    }

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

    std::ofstream out("force_nodes.csv");
    out << "x,y,z\n";

    for (auto& node : nodes) 
    {
        if (node->status == NodeStatus::CUT)
        {
            continue;
        }

        Point p = node->position;

        CGAL::Bounded_side side = inside_test(p);

        /*
         * Trying a new way to pick force nodes where if you are inside of this then you are a force node because cut nodes are already cut out
        if (side == CGAL::ON_BOUNDED_SIDE && side != CGAL::ON_BOUNDARY)
        {
            node->status = NodeStatus::CUT_NEXT;
        } 

        else if (side == CGAL::ON_BOUNDARY) 
        {
            if (side == CGAL::ON_BOUNDARY)
            {
                node->status = NodeStatus::FORCE_NODE;
                out << node->position.x() << ',' << node->position.y() << ',' << node->position.z() << '\n';

            }
        }
        */

        if (side == CGAL::ON_BOUNDED_SIDE || side == CGAL::ON_BOUNDARY)
        {
                node->status = NodeStatus::FORCE_NODE;
                out << node->position.x() << ',' << node->position.y() << ',' << node->position.z() << '\n';
        }

    }
}



/*
 * Takes in a tool path and returns a vector of all of the elements that will be removed from this cut
*/
std::vector<Tr::Cell_handle> MeshHandler::set_simulated_tool_path(const Surface_mesh& tool_path_mesh)
{
    Tree tree(faces(tool_path_mesh).first, faces(tool_path_mesh).second, tool_path_mesh);
    tree.accelerate_distance_queries();

    CGAL::Side_of_triangle_mesh<Surface_mesh, Kernel> inside_test(tree);

    std::vector<Tr::Cell_handle> removed_elems;
    const auto& tr = volume_mesh.triangulation();
    for (auto cit = tr.finite_cells_begin(); cit != tr.finite_cells_end(); ++cit)
    {
        Point v1 = cit->vertex(0)->point().point();
        Point v2 = cit->vertex(1)->point().point();
        Point v3 = cit->vertex(2)->point().point();
        Point v4 = cit->vertex(3)->point().point();
        
        Point centroid ((v1.x() + v2.x() + v3.x() + v4.x()) / 4,
                        (v1.y() + v2.y() + v3.y() + v4.y()) / 4,
                        (v1.z() + v2.z() + v3.z() + v4.z()) / 4);

;

        CGAL::Bounded_side side = inside_test(centroid);

        /*
         * Trying a new way to pick force nodes where if you are inside of this then you are a force node because cut nodes are already cut out
        if (side == CGAL::ON_BOUNDED_SIDE && side != CGAL::ON_BOUNDARY)
        {
            node->status = NodeStatus::CUT_NEXT;
        } 

        else if (side == CGAL::ON_BOUNDARY) 
        {
            if (side == CGAL::ON_BOUNDARY)
            {
                node->status = NodeStatus::FORCE_NODE;
                out << node->position.x() << ',' << node->position.y() << ',' << node->position.z() << '\n';

            }
        }
        */

        if (side == CGAL::ON_BOUNDED_SIDE || side == CGAL::ON_BOUNDARY)
        {
            removed_elems.push_back(cit);
        }

    }

    return removed_elems;
}



/*
 * Custom function to get the intersection volume between a volumetric mesh and a surface mesh because CGAL doesn't provide one
 * Works by testing for every element in the volumetric mesh to see if it's in the surface mesh and sum up all the elements that are
*/
static double intersection_volume(const C3t3& c3,
                           Surface_mesh& shell)
{
    /*
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

    double shell_volume = CGAL::Polygon_mesh_processing::volume(shell);

    std::cout << "Shell volume: " << shell_volume << std::endl;
    std::cout << "intersection volume: " << intersection_volume << std::endl;
    
    return CGAL::to_double(intersection_volume) / 
    (shell_volume + CGAL::to_double(total_volume) - CGAL::to_double(intersection_volume)); // intersection / Union
    */
    
    namespace PMP = CGAL::Polygon_mesh_processing;

    Kernel::FT intersection_volume = 0.0;
    Kernel::FT total_volume = 0;
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

        total_volume += PMP::volume(tet_mesh);
        PMP::clip(clipped, shell_copy, PMP::parameters::clip_volume(true));

        if (!clipped.is_empty())
            intersection_volume += PMP::volume(clipped);
    }

    double shell_volume = CGAL::Polygon_mesh_processing::volume(shell);

    std::cout << "Shell volume: " << shell_volume << std::endl;
    std::cout << "intersection volume: " << intersection_volume << std::endl;
    
    return CGAL::to_double(intersection_volume) / 
    (shell_volume + CGAL::to_double(total_volume) - CGAL::to_double(intersection_volume)); // intersection / Union


    /*
    Surface_mesh sm;
    CGAL::facets_in_complex_3_to_triangle_mesh(c3, sm);

    // 1) Repair topological issues
    CGAL::Polygon_mesh_processing::stitch_borders(sm);

    // 2) Remove degenerate / collapsed triangles
    CGAL::Polygon_mesh_processing::remove_degenerate_faces(sm);

    // 3) Remove isolated vertices
    CGAL::Polygon_mesh_processing::remove_isolated_vertices(sm);

    // 4) Remove self-intersections
    std::vector<std::pair<CGAL::SM_Face_index, CGAL::SM_Face_index>> self_ints;
    CGAL::Polygon_mesh_processing::self_intersections(sm, std::back_inserter(self_ints));
    for (auto& f : self_ints)
    {
        sm.remove_face(f.first);
    }

    CGAL::Polygon_mesh_processing::triangulate_faces(sm);
    CGAL::Polygon_mesh_processing::isotropic_remeshing(
        faces(sm),
        0.3,
        sm
    );

    // 5) Orient consistently
    CGAL::Polygon_mesh_processing::orient(sm);

    // 6) Check watertight
    if (!CGAL::is_closed(sm))
    {
        std::cerr << "WARNING: mesh is not closed — corefinement will fail!\n";
    }


    Surface_mesh intersection_mesh;
    std::cout << "BEFORE interseciont computation" << std::endl;
    CGAL::Polygon_mesh_processing::corefine_and_compute_intersection(sm, shell, intersection_mesh);
    std::cout << "AFTER interseciont computation" << std::endl;

    return {CGAL::to_double(CGAL::Polygon_mesh_processing::volume(intersection_mesh)), CGAL::to_double(CGAL::Polygon_mesh_processing::volume(sm))};
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
    //double goal_volume = context->goal_volume;

    double pressure = x[0];

    FEA fea(*handler);

    C3t3 deformed_mesh = fea.apply_trench_force(pressure);

    double intersection_over_union = intersection_volume(deformed_mesh, *goal_mesh);
    std::cout << "[OBJ] pressure=" << pressure << '\n'
    << "Intersection over Union=" << intersection_over_union
    << "\n";

    export_surface_mesh("debug_mesh.stl", &deformed_mesh);
    return intersection_over_union; // NLopt will maximize this

}


/*
 * This function takes in a tool path mesh of the cut and the deformed mesh after the cut and calculates the stress
 * It stores the stress values in this class so it doesn't return anything
*/
void MeshHandler::set_stress_in_cut(const std::filesystem::path& tool_path_stl, const std::filesystem::path& real_mesh_stl)
{
    export_surface_mesh("starting_mesh.stl", &volume_mesh);
    set_tool_path(tool_path_stl); // sets all values in nodes vector to what we need to simulate the cut
    
    std::cout << "Summary after tool path" << std::endl;
    print_summary();

    Surface_mesh real_mesh;

    if (!CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(real_mesh_stl.string(), real_mesh)) 
    {
        throw std::runtime_error("Failed to load mesh: " + real_mesh_stl.string());
    }

    OptimizationContext context{this, &real_mesh};

    nlopt::opt opt(nlopt::LN_COBYLA, 1);  // 1 parameter: pressure
    opt.set_max_objective(objective_func, &context);
    //opt.set_xtol_rel(1e-3);
    opt.set_ftol_rel(1e-4);
    opt.set_maxeval(40);

    std::vector<double> x = {1}; // NOTE: x[0] is pressure. We just use it as a vector because that's what NLOPT likes.
    double max_volume;
    opt.optimize(x, max_volume);


    std::cout << "[Optimization] Best pressure = "
              << x[0] << " MPa (volume = " << max_volume << ")\n";



    FEA fea(*this);

    fea.get_cut_stress_field(x[0]);

    // Now that we have set the stress in the cut we need to make the mesh match exactly to the real mesh
    // After matching the pressure our mesh will almost certainly not be an identical match and require us to match it
    // To do this we map all of the required deformations to map every exterior node in our working mesh onto the exterior of the 
    // real mesh scan. We can then run an FEA analysis which will get us the interior mapping that has the lowest elastic strain 
    // left inside of the material. This doesn't give us a perfect mapping but it's about as good as we can do.

    std::vector<std::pair<size_t, MyVector>> displacements; // this stores all of the exterior nodes we need to pass to our FEA function
    Tree real_mesh_tree(faces(real_mesh).begin(), faces(real_mesh).end(), real_mesh);
    real_mesh_tree.accelerate_distance_queries();

    for (size_t i = 0; i < nodes.size(); ++i)
    {
        std::shared_ptr<Node> node = nodes[i];

        if (node->status == NodeStatus::CUT)
        {
            continue;
        }

        if (node->status == NodeStatus::CUT_NEXT)
        {
            node->status = NodeStatus::CUT;
            continue;
        }

        if (node->status == NodeStatus::FORCE_NODE)
        {
            node->status = NodeStatus::KNOWN_WALL;
            node->pressure = x[0];
            node->is_exterior = true; // If we just applied force to you you're now on the exterior
        }

        if (node->is_exterior)
        {
            Point point = real_mesh_tree.closest_point(node->position);
            MyVector displacement = {point.x() - node->position.x(),
                                    point.y() - node->position.y(),
                                    point.z() - node->position.z()};
                                    
            displacements.emplace_back(i, displacement);
        }
    }


    std::cout << "Calling displace_enterior, please don't fail" << std::endl;
    volume_mesh = fea.displace_enterior(displacements);

    export_surface_mesh("displaced_exterior_mesh.stl", &volume_mesh);

}


C3t3 MeshHandler::simulate_cut(const Surface_mesh& tool_path)
{
    std::vector<Tr::Cell_handle> cut_elems = set_simulated_tool_path(tool_path);

    FEA fea(*this);
    return fea.simulate_cut(cut_elems);
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
    size_t cut_next_count = 0;
    size_t cut_count = 0;
    size_t known_wall_count = 0;

    for (const auto& node : nodes) 
    {
        if (node->is_boundary) boundary_count++;

        switch (node->status) 
        {
            case NodeStatus::UN_CUT:    ++uncut_count; break;
            case NodeStatus::CUT_NEXT:  ++cut_next_count; break;
            case NodeStatus::CUT:       ++cut_count; break;
            case NodeStatus::FORCE_NODE: ++force_count; break;
            case NodeStatus::KNOWN_WALL: ++known_wall_count; break;
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
    std::cout << "   CUT_NEXT:      " << cut_next_count << " (" << pct(cut_next_count) << "%)\n";
    std::cout << "   CUT:           " << cut_count     << " (" << pct(cut_count)     << "%)\n";
    std::cout << "   KNOWN_WALL:    " << known_wall_count << " (" << pct(known_wall_count) << "%)\n";
    std::cout << "==================================\n\n";
}


/*
 * Outputs the current volumetric mesh stored in the class as an stl file
 * This is very helpful for testing
*/
void export_surface_mesh(const std::filesystem::path& out_path, const C3t3 *input_mesh)
{
    Surface_mesh surface_mesh;
    CGAL::facets_in_complex_3_to_triangle_mesh(*input_mesh, surface_mesh);

    if (!CGAL::IO::write_polygon_mesh(out_path.string(), surface_mesh, CGAL::parameters::stream_precision(17))) {
        throw std::runtime_error("Failed to write surface mesh to " + out_path.string());
    }

    std::cout << "[MeshHandler] Exported surface mesh to " << out_path << std::endl;
}
