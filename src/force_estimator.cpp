#include "force_estimator.hpp"

struct RoundedPt 
{
    long long x, y, z;

    // C++17-compatible equality:
    bool operator==(const RoundedPt& other) const noexcept
    {
        return x == other.x && y == other.y && z == other.z;
    }
};

struct RoundedHash 
{
    std::size_t operator()(const RoundedPt& p) const noexcept
    {
        // simple XOR-combine hash
        std::size_t h1 = CGAL::hash_value(p.x); std::size_t h2 = CGAL::hash_value(p.y);
        std::size_t h3 = CGAL::hash_value(p.z);
        return h1 ^ (h2 << 1) ^ (h3 << 2);
    }
};


constexpr double EPS = 1e-8;               // rounding grid
auto key = [](const Point& p)
{
    return RoundedPt{ std::llround(p.x()/EPS),
                      std::llround(p.y()/EPS),
                      std::llround(p.z()/EPS) };
};

ToolPath::ToolPath(const std::filesystem::path& stl_file)
{
	assert(CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(stl_file.string(), tool_path_mesh));
}


std::vector<Node> ToolPath::create_nodes(C3t3& mesh_3d)
{
	Surface_mesh mesh;
	CGAL::facets_in_complex_3_to_triangle_mesh(mesh_3d, mesh);
	std::unordered_map<Vertex_descriptor, Vector> vertex_normals;
	boost::associative_property_map<std::unordered_map<Vertex_descriptor, Vector>> vertex_normals_map(vertex_normals);

	if (!CGAL::Polygon_mesh_processing::is_outward_oriented(mesh)) 
	    CGAL::Polygon_mesh_processing::reverse_face_orientations(mesh);  

	CGAL::Polygon_mesh_processing::compute_vertex_normals(
						mesh, 
						vertex_normals_map,
						CGAL::Polygon_mesh_processing::parameters::vertex_point_map(get(CGAL::vertex_point, mesh))
						);
	// Build the AABB tree on the bounding mesh
	CGAL::AABB_tree<CGAL::AABB_traits_3<Kernel, CGAL::AABB_face_graph_triangle_primitive<Surface_mesh>>> 
	AABB_tree(faces(tool_path_mesh).begin(), faces(tool_path_mesh).end(), tool_path_mesh);

	AABB_tree.accelerate_distance_queries();

	std::vector<Node> result;
	std::vector<Vertex_descriptor> bottom_vertices;
	std::unordered_map<RoundedPt, Node, RoundedHash> point_to_node;
	std::cout << "# of vertices in mesh " << mesh.vertices().size() << '\n';


	double tolerance = 0.001;
	double tolerance_squared = tolerance * tolerance;
	for (const Vertex_descriptor& vertex : vertices(mesh))
	{
		Point point = mesh.point(vertex);
		Vector normal_vector = vertex_normals[vertex];

		double distance_squared = AABB_tree.squared_distance(point);

		if (distance_squared < tolerance_squared)
		{
			if (std::abs(normal_vector.z()) <= std::max(std::abs(normal_vector.x()), std::abs(normal_vector.y())))
			{
				double length = std::sqrt(normal_vector.x()*normal_vector.x() + normal_vector.y()*normal_vector.y());
				UnitVector unit_vector = {normal_vector.x() / length, normal_vector.y() / length};
				result.emplace_back(vertex, point, unit_vector);
				point_to_node[key(point)] = result.back();
			}
			else	
			{
				bottom_vertices.push_back(vertex);
			}
		}

	}

	std::cout << "# of kd_tree nodes: " << result.size() << std::endl;

	CGAL::Kd_tree<CGAL::Search_traits_3<Kernel>> kd_tree;
	for (const auto& node : result)
	{
		kd_tree.insert(node.position);
	}

	for (const Vertex_descriptor& vertex : bottom_vertices)
	{
		Point point = mesh.point(vertex);
		CGAL::K_neighbor_search<CGAL::Search_traits_3<Kernel>> search(kd_tree, point, 1);
		Point closest_point = search.begin()->first;
		assert(point_to_node.find(key(closest_point)) != point_to_node.end()); result.emplace_back(vertex, point, point_to_node[key(closest_point)].normal_vector);
	}

	return result;
}

static double tet_quality(const Point& a, const Point& b, const Point& c, const Point& d) 
{
    // volume
    double V = std::abs(CGAL::to_double(CGAL::volume(a,b,c,d)));
    // average edge length
    double e2 = 0; int m=0;
    Point P[4] = {a,b,c,d};
    for(int i=0;i<4;i++) for(int j=i+1;j<4;j++){ 
	e2 += CGAL::to_double(CGAL::squared_distance(P[i],P[j])); ++m; 
    }
    double e_avg = std::sqrt(e2/m);
    return V / (e_avg*e_avg*e_avg); // dimensionless; ~0 for slivers
}


ResidualCalculator::ResidualCalculator(const std::filesystem::path& object_stl, const std::filesystem::path& tool_path_stl, 
		    const std::filesystem::path& boundary_stl, const std::filesystem::path& goal_stl)
{
	Surface_mesh boundary_mesh;
	Surface_mesh pre_deformation_surface_mesh;
	assert(CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(object_stl.string(), pre_deformation_surface_mesh));
	assert(CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(goal_stl.string(), goal_mesh)); 
	assert(CGAL::Polygon_mesh_processing::IO::read_polygon_mesh(boundary_stl.string(), boundary_mesh));

	Domain domain(pre_deformation_surface_mesh);
	//domain.detect_features();

	Criteria crit(CGAL::parameters::edge_size(5)
			      .facet_angle(5)
			      .facet_distance(1)
			      .cell_radius_edge_ratio(3.0)
			      .cell_size(10));

	std::cout << "# of faces: " << pre_deformation_surface_mesh.num_faces() << std::endl;
	mesh_c3t3 = CGAL::make_mesh_3<C3t3>(domain, crit);


	std::cout << "# of tetrahedrals: " << mesh_c3t3.number_of_cells() << std::endl;

	ToolPath tool_path(tool_path_stl);
	force_nodes = tool_path.create_nodes(mesh_c3t3);

	for (auto vit = mesh_c3t3.triangulation().finite_vertices_begin(); vit != mesh_c3t3.triangulation().finite_vertices_end(); ++vit)
	{
		point_to_vertex_handle[vit->point().point()] = vit;
	}


	std::cout << "# of force nodes " << force_nodes.size() << '\n';

	size_t id = 1;
	// Access the triangulation from the C3t3 object
	const auto& triangulation = mesh_c3t3.triangulation();

	for (auto vit = triangulation.finite_vertices_begin(); vit != triangulation.finite_vertices_end(); ++vit) 
	{
		node_id_map[vit] = id;
		id_node_map[id] = vit;
		id_point_map[id] = vit->point().point();
		id++;
	}

	auto mid_id = id_point_map.size(); // start allocating after corner nodes
	auto get_mid = [&](size_t a, size_t b)->size_t 
	{

	    if (a>b) 
	    {
		std::swap(a,b);
	    }
	    Edge E{a,b};
	    auto it = edge_mid_map.find(E);

	    if (it != edge_mid_map.end()) 
	    {
		    return it->second;
	    }
	    ++mid_id;
	    // midpoint
	    Point p = Point( (id_point_map[a].x()+id_point_map[b].x())*0.5,
			     (id_point_map[a].y()+id_point_map[b].y())*0.5,
			     (id_point_map[a].z()+id_point_map[b].z())*0.5 );
	    id_point_map[mid_id] = p; // grow vector so we can later write *NODE
	    edge_mid_map[E] = mid_id;
	    return mid_id;
	};


	// ---- Precompute C3D10 connectivity for all elements ----
	elem10.reserve(mesh_c3t3.triangulation().number_of_finite_cells());

	// Build the AABB tree on the boundary mesh
	CGAL::AABB_tree<CGAL::AABB_traits_3<Kernel, CGAL::AABB_face_graph_triangle_primitive<Surface_mesh>>> 
	tree(faces(boundary_mesh).begin(), faces(boundary_mesh).end(), boundary_mesh);
	tree.accelerate_distance_queries(); // Construct side tester
	CGAL::Side_of_triangle_mesh<Surface_mesh, Kernel> inside_test(tree);
	for (auto vit = mesh_c3t3.triangulation().finite_vertices_begin(); vit != mesh_c3t3.triangulation().finite_vertices_end(); ++vit)
	{
		Point point = vit->point().point();
		if (inside_test(point) == CGAL::ON_BOUNDED_SIDE)
		{
			boundary_nodes.push_back(vit);
		}
	}

	std::unordered_set<Tr::Vertex_handle> c3t3_force_nodes;

	for (const auto& node : force_nodes)
	{
		Tr::Vertex_handle handle = point_to_vertex_handle[node.position];
		c3t3_force_nodes.insert(handle);
	}


	std::size_t eid = 1;
	for (auto cit = mesh_c3t3.triangulation().finite_cells_begin();
	     cit != mesh_c3t3.triangulation().finite_cells_end(); ++cit, ++eid)
	{


	    // store connectivity (for later *ELEMENT block)
	    std::array<std::size_t,4> conn;
	    for(int v=0; v<4; ++v) 
	    {
		 conn[v] = node_id_map[cit->vertex(v)];
	    }

	    if (tet_quality(id_point_map[conn[0]], id_point_map[conn[1]], id_point_map[conn[2]], id_point_map[conn[3]]) < 0.05) continue;

	    // examine each of the 4 faces
	    for(int opp=0; opp<4; ++opp)
	    {
		// build the 3 vertices that make up this face
		Tr::Vertex_handle v0 = cit->vertex((opp+1)&3);
		Tr::Vertex_handle v1 = cit->vertex((opp+2)&3);
		Tr::Vertex_handle v2 = cit->vertex((opp+3)&3);

		// if *all three* belong to the trench-node set, this face is on wall
		if (c3t3_force_nodes.count(v0) &&
		    c3t3_force_nodes.count(v1) &&
		    c3t3_force_nodes.count(v2))
		{
		    trench_faces.emplace_back(eid, opp+1);   // opp==0 → S1, etc.
		}
	    }


	    //now build C3D10 faces for calculix use
	    size_t n1 = node_id_map[cit->vertex(0)];
	    size_t n2 = node_id_map[cit->vertex(1)];
	    size_t n3 = node_id_map[cit->vertex(2)];
	    size_t n4 = node_id_map[cit->vertex(3)];

	    size_t n12 = get_mid(n1,n2);
	    size_t n23 = get_mid(n2,n3);
	    size_t n31 = get_mid(n3,n1);
	    size_t n14 = get_mid(n1,n4);
	    size_t n24 = get_mid(n2,n4);
	    size_t n34 = get_mid(n3,n4);

	    //check for orientation
	    auto p0 = cit->vertex(0)->point().point();
	    auto p1 = cit->vertex(1)->point().point();
	    auto p2 = cit->vertex(2)->point().point();
	    auto p3 = cit->vertex(3)->point().point();
	    if (CGAL::orientation(p0,p1,p2,p3) != CGAL::POSITIVE) 
	    {
		// swap to fix orientation (e.g., swap v1<->v2)
		std::swap(n2, n3);
		std::swap(p1, p2);
	    }

	    elem10.push_back({ n1,n2,n3,n4, n12,n23,n31,n14,n24,n34 });
	    
	}

	std::cout << "# of boundary nodes " << boundary_nodes.size() << '\n';
}


template<class K>
CGAL::Aff_transformation_3<K>
make_rigid_transform(double t_x,double t_y,double t_z,
                     double r_x,double r_y,double r_z)
{
    using FT = typename K::FT;
    FT cx = std::cos(r_x), sx = std::sin(r_x);
    FT cy = std::cos(r_y), sy = std::sin(r_y);
    FT cz = std::cos(r_z), sz = std::sin(r_z);

    CGAL::Aff_transformation_3<K> Rx(  1,0,0,0,
                                        0,cx,-sx,0,
                                        0,sx, cx,0, 1 );

    CGAL::Aff_transformation_3<K> Ry(  cy,0,sy,0,
                                        0, 1, 0,0,
                                       -sy,0,cy,0, 1 );

    CGAL::Aff_transformation_3<K> Rz(  cz,-sz,0,0,
                                        sz, cz,0,0,
                                         0,  0,1,0, 1 );

    CGAL::Aff_transformation_3<K> T(CGAL::TRANSLATION,
                                    typename K::Vector_3(t_x,t_y,t_z));

    return  T * Rz * Ry * Rx;   // translate ∘ Rz ∘ Ry ∘ Rx
}

// ── generic helper that works for any Triangulation_3 or C3t3 ────────────────
template<class C3, class Aff>
void transform_volumetric_mesh(C3& c3, const Aff& M)
{
    auto& tr = c3.triangulation();

    for(auto vit = tr.finite_vertices_begin();
        vit != tr.finite_vertices_end(); ++vit)
    {
        const CGAL::Weighted_point_3<Kernel>& wp  = vit->point();              // (p, w)
        const Point&     p   = wp.point();                // strip weight
        const Point      p2  = M(p);                      // affine image

        vit->point() = CGAL::Weighted_point_3<Kernel>(p2, wp.weight());        // keep weight
    }
}

// ── your new API: works for C3t3; call just like the old rotate_mesh() ───────
template<class C3t3>
void rotate_mesh(C3t3& mesh,
                 double t_x,double t_y,double t_z,
                 double r_x,double r_y,double r_z)
{
    using Kernel = typename C3t3::Triangulation::Geom_traits;
    auto M = make_rigid_transform<Kernel>(t_x,t_y,t_z,r_x,r_y,r_z);
    transform_volumetric_mesh(mesh, M);
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


double ResidualCalculator::obj(double pressure, double t_x, double t_y, double t_z,
			       double r_x, double r_y, double r_z)
{
	std::cout << "pressure: " << pressure << std::endl;
	std::cout << "t_x: " << t_x << std::endl;
	std::cout << "t_y: " << t_y << std::endl;
	std::cout << "t_z: " << t_z << std::endl;
	std::cout << "r_x: " << r_x << std::endl;
	std::cout << "r_y : " << r_y << std::endl;
	std::cout << "r_z: " << r_z<< std::endl;

	const auto& triangulation = mesh_c3t3.triangulation();

	std::string calculix_input_file {"calculix"};

	std::ofstream inp(calculix_input_file + ".inp");

	inp << "*NODE\n";
	for (auto vit = triangulation.finite_vertices_begin(); vit != triangulation.finite_vertices_end(); ++vit) 
	{
		size_t id = node_id_map[vit];
		const Point &point = vit->point().point();
		inp << id << ", " << point.x() << ", " << point.y() << ", " << point.z() << '\n';
	}

	const size_t corner_count = node_id_map.size(); 
	for (size_t id = corner_count+1; id <= id_point_map.size(); ++id) 
	{
	    const Point& p = id_point_map[id];
	    inp << id << ", " << p.x() << ", " << p.y() << ", " << p.z() << '\n';
	}


	inp << "*NSET, NSET=ALL, GENERATE\n";
	inp << 1 << ", " << id_point_map.size() << ", " << 1 << '\n';
	inp << "*ELEMENT, TYPE=C3D10, ELSET=TRENCH_ELEMS\n";

	for(size_t i=0; i < elem10.size() ; ++i)
	{
	    const auto& c = elem10[i];
	    inp << i + 1;
	    for (int k=0; k<10; ++k) 
	    {
		inp << ", " << c.v[k];
	    }
	    inp << '\n';
	}

	inp << "*SURFACE, TYPE=ELEMENT, NAME=TRENCH_FACES\n";
	std::unordered_set<size_t> layer1_ids;
	for(const auto& ef : trench_faces)
	{
	    inp << ef.first << ", S" << ef.second << '\n';
	    layer1_ids.insert(ef.first);
	}

	inp << "*ELSET, ELSET=LAYER1\n";
	std::vector<size_t> ids(layer1_ids.begin(), layer1_ids.end());
	std::sort(ids.begin(), ids.end());
	const int per_line = 16;
	for (size_t i = 0; i < ids.size(); i += per_line) 
	{
	    size_t end = std::min(i + per_line, ids.size());
	    for (size_t j = i; j < end; ++j) 
	    {
		inp << ids[j];
		if (j + 1 < end) inp << ", ";
	    }
	    inp << "\n";
	}
	
	inp << "*MATERIAL, NAME=ALUMINUM\n";
	inp << "*ELASTIC\n";
	inp << "69000, 0.33\n";  // MPa units

	inp << "*SOLID SECTION, ELSET=TRENCH_ELEMS, MATERIAL=ALUMINUM\n";

	inp << "*STEP\n";
	inp << "*STATIC\n";

	inp << "*BOUNDARY\n";

	for (const auto &vit : boundary_nodes) 
	{
		size_t id = node_id_map[vit];
		inp << id << ", 1, 3, 0.0\n";
	}

	// Inside the step:
	inp << "*DLOAD\n";
	inp << "TRENCH_FACES, P, " << pressure << "\n";  // p is your scalar decision variable (MPa)

	inp << "*NODE PRINT, NSET=ALL\n";
	inp << "U\n";
	inp << "*EL PRINT, ELSET=LAYER1, FREQUENCY=0\n";
	inp << "  S\n";
	inp << "*END STEP\n";

	inp.close();

	std::string run_calculix {"ccx " + calculix_input_file};
	run_calculix += " > /dev/null 2>&1";
	int result = std::system(run_calculix.c_str());
	assert(result == 0);

	std::ifstream output("calculix.dat");
	if (!output.is_open()) 
	{
	    throw std::runtime_error("Unable to open calculix.dat");
	}

	bool in_displacement = false;
	std::string line;
	size_t node_id;
	double delta_x;
	double delta_y;
	double delta_z;
	std::unordered_map<size_t, Vector> id_to_displacement;
	while (getline(output, line))
	{
		if (line.find("displacement") != std::string::npos)
		{
			in_displacement = true;
			continue;
		}

		if (in_displacement && line.empty())
		{
			continue;
		}

		if (in_displacement)
		{
			std::size_t id; double dx,dy,dz;
			std::istringstream(line) >> id >> dx >> dy >> dz;
			id_to_displacement[id] = Vector(dx,dy,dz);
		}
	}

	C3t3 test_mesh = mesh_c3t3;                // deep copy

	//------------------------------------------------------
	// Build id -> vertex-handle map *for the copy*
	// (iterate in the *same* order you used when writing *NODE)
	//------------------------------------------------------
	std::unordered_map<std::size_t, Tr::Vertex_handle> id_to_vh_copy;
	std::size_t id = 1;
	for(auto vit = test_mesh.triangulation().finite_vertices_begin();
	    vit != test_mesh.triangulation().finite_vertices_end(); ++vit, ++id)
	{
	    id_to_vh_copy[id] = vit;
	}

	//------------------------------------------------------
	// Apply the displacements
	//------------------------------------------------------
	for(const auto& kv : id_to_displacement) 
	{
	    std::size_t id         = kv.first;
	    const Vector& disp     = kv.second;

	    auto hit = id_to_vh_copy.find(id);
	    if(hit == id_to_vh_copy.end()) continue;      // safety

	    Tr::Vertex_handle vh = hit->second;
	    Point  p  = vh->point().point();
	    vh->point() = Tr::Point(p + disp);
	}
	rotate_mesh(test_mesh, t_x, t_y, t_z, r_x, r_y, r_z);
		

	double volume = intersection_volume(test_mesh, goal_mesh);
	std::cout << "Volume: " << volume << std::endl;
	return volume;
		
}


static double my_objective(const std::vector<double> &x, std::vector<double> &grad, void *data) 
{
	ResidualCalculator *rc = static_cast<ResidualCalculator*>(data);
	return rc->obj(x[0], x[1], x[2], x[3], x[4], x[5], x[6]);
}


std::unordered_map<Tr::Vertex_handle, StressTensor> ResidualCalculator::get_stress_tensor()
{
	std::unordered_map<Tr::Vertex_handle, StressTensor> handle_to_tensor;
	std::ifstream frd("calculix.frd");
	if (!frd) throw std::runtime_error("cannot open calculix.frd");

	std::string line;
	bool in_frame = false;

	while (std::getline(frd, line))
	{
		// toggle on/off inside the single result frame
		if (line.find("3C") == 0) 
		{
			std::cout << "in frame now" << '\n';
			in_frame = !in_frame; 
			continue; 
		}

		if (!in_frame) 
		{
			continue;
		}

		// node records start with " 1"
		if (line.size() < 2 || line[0]!=' ' || line[1]!='1') 
		{
			continue;
		}

		// expected layout:
		//  1 <nodeID> 1  S  <Sxx> <Syy> <Szz> <Sxy> <Sxz> <Syz>
		std::istringstream iss(line);
		int rec,node,one; char key;
		std::array<double,6> s{};
		iss >> rec >> node >> one >> key
		    >> s[0] >> s[1] >> s[2] >> s[3] >> s[4] >> s[5];

		auto vit = id_node_map.find(node);
		if (vit != id_node_map.end())
		{
		    handle_to_tensor[vit->second] = StressTensor{s[0], s[1], s[2]};   // store by Vertex_handle
		}
	}

	return handle_to_tensor;
}


std::unordered_map<Tr::Vertex_handle, StressTensor> ResidualCalculator::stress_estimator()
{
	nlopt::opt optimizer(nlopt::LN_COBYLA, 7);
	optimizer.set_max_objective(my_objective, this);
	optimizer.set_lower_bounds({-5000, -5, -5, -5, -M_PI, -M_PI, -M_PI});
	optimizer.set_upper_bounds({5000,  5, 5,  5,  M_PI,  M_PI,  M_PI});

	optimizer.set_initial_step({10, 0.5, 0.5, 0.5, 0.1, 0.1, 0.1});

	optimizer.set_xtol_rel(1e-3);
	std::vector<double> x = {0, 0, 0, 0, 0, 0, 0};
	double max_intersection;
	optimizer.optimize(x, max_intersection);
	std::unordered_map<Tr::Vertex_handle, StressTensor> handle_to_tensor = get_stress_tensor();
	return handle_to_tensor;
}
