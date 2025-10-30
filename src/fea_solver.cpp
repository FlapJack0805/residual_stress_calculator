#include "fea_solver.hpp"

static bool face_compare(const FaceRef& a, const FaceRef& b)
{
    if (a.elem_id < b.elem_id) return true;
    if (a.elem_id > b.elem_id) return false;
    return a.face_id < b.face_id;
}

FEA::FEA(MeshHandler& mesh_in) : mesh(mesh_in) {}


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


/*
 * This takes in a pressure and applys that pressure to the trench walls
 * It outputs a new C3t3 mesh with the deformation caused by said pressure
*/
C3t3 FEA::apply_trench_force(double pressure)
{
    const auto& nodes = mesh.get_nodes();
    const auto& c3t3 = mesh.get_c3t3();
    const auto& tr = c3t3.triangulation();

    std::string BASE_FILE_NAME = "calculix";

    std::ofstream inp(BASE_FILE_NAME + ".inp");
    if (!inp.is_open())
        throw std::runtime_error("Failed to open CalculiX input file");

    // -------------------------------
    // *NODE block
    // -------------------------------
    inp << "*NODE\n";
    std::map<Tr::Vertex_handle, size_t> node_id_map;
    size_t id = 1;

    for (auto vit = tr.finite_vertices_begin(); vit != tr.finite_vertices_end(); ++vit, ++id) {
        const Point p = vit->point().point();
        node_id_map.insert({vit, id});
        inp << id << ", " << p.x() << ", " << p.y() << ", " << p.z() << "\n";
    }

    inp << "*NSET, NSET=ALL, GENERATE\n1, " << id-1 << ", 1\n";

    // -------------------------------
    // *ELEMENT block (C3D10)
    // -------------------------------
    //inp << "*ELEMENT, TYPE=C3D10, ELSET=TRENCH_ELEMS\n";
    //will switch to C3D10 later because it's better but C3D4 is easier to get this working first
    inp << "*ELEMENT, TYPE=C3D4, ELSET=TRENCH_ELEMS\n";
    size_t elem_id = 1;

    for (auto cit = tr.finite_cells_begin(); cit != tr.finite_cells_end(); ++cit, ++elem_id) {

        if (tet_quality(cit->vertex(0)->point().point(), cit->vertex(1)->point().point(), cit->vertex(2)->point().point(), cit->vertex(3)->point().point()) < 0.05)
        {
            continue;
        }

        inp << elem_id;
        for (int i = 0; i < 4; ++i)
            inp << ", " << node_id_map.at(cit->vertex(i));
        inp << "\n";
    }

    // -------------------------------
    // MATERIAL and SECTION
    // -------------------------------
    inp << "*MATERIAL, NAME=ALUMINUM\n";
    inp << "*ELASTIC\n";
    inp << "69000, 0.33\n"; // MPa
    inp << "*SOLID SECTION, ELSET=TRENCH_ELEMS, MATERIAL=ALUMINUM\n";

    // -------------------------------
    // BOUNDARY CONDITIONS
    // -------------------------------
    inp << "*BOUNDARY\n";
    for (const auto& node : nodes) {
        if (node->is_boundary) {
            inp << node_id_map.at(node->handle) << ", 1, 3, 0.0\n";
        }
    }

    // -------------------------------
    // SURFACE SET for trench wall
    // -------------------------------
    // Collect trench faces
    std::vector<FaceRef> trench_faces;

    size_t elem_counter = 1;
    for (auto cit = tr.finite_cells_begin(); cit != tr.finite_cells_end(); ++cit, ++elem_counter)
    {
        for (int face = 0; face < 4; ++face)
        {
            int v0 = (face + 1) % 4;
            int v1 = (face + 2) % 4;
            int v2 = (face + 3) % 4;

            auto vh0 = cit->vertex(v0);
            auto vh1 = cit->vertex(v1);
            auto vh2 = cit->vertex(v2);

            const auto& n0 = nodes[node_id_map.at(vh0) - 1];
            const auto& n1 = nodes[node_id_map.at(vh1) - 1];
            const auto& n2 = nodes[node_id_map.at(vh2) - 1];

            // trying or instead of and for a bit
            if (n0->status == NodeStatus::FORCE_NODE ||
                n1->status == NodeStatus::FORCE_NODE ||
                n2->status == NodeStatus::FORCE_NODE)
            {
                Vector N = CGAL::cross_product(n1->position - n0->position,
                                               n2->position - n0->position);

                double ax = std::abs(N.x());
                double ay = std::abs(N.y());
                double az = std::abs(N.z());

                // Skip horizontal faces (floor)
                //if (az > ax && az > ay) continue;
                //std::cout << "adding load to face" << std::endl;

                trench_faces.push_back({elem_counter, face + 1});
            }
        }
    }

    // ---- Sort faces for consistent CalculiX output ----
    std::sort(trench_faces.begin(), trench_faces.end(), face_compare);

    // ---- Write SURFACE section ----
    inp << "*SURFACE, TYPE=ELEMENT, NAME=TRENCH_FACES\n";
    for (auto& f : trench_faces)
        inp << f.elem_id << ", S" << f.face_id << "\n";

    // ---- LAYER1 ELSET (unique element IDs only) ----
    inp << "*ELSET, ELSET=LAYER1\n";

    std::unordered_set<size_t> unique_elems;
    for (auto& f : trench_faces) unique_elems.insert(f.elem_id);

    std::vector<size_t> elem_ids(unique_elems.begin(), unique_elems.end());
    std::sort(elem_ids.begin(), elem_ids.end());

    const int per_line = 16;
    for (size_t i = 0; i < elem_ids.size(); i += per_line)
    {
        size_t end = std::min(i + per_line, elem_ids.size());
        for (size_t j = i; j < end; ++j)
        {
            inp << elem_ids[j];
            if (j + 1 < end) inp << ", ";
        }
        inp << "\n";
    }

    // -------------------------------
    // LOAD STEP (apply pressure)
    // -------------------------------
    inp << "*STEP\n";
    inp << "*STATIC\n";
    inp << "*DLOAD\n";
    inp << "TRENCH_FACES, P, " << pressure << "\n";
    inp << "*NODE PRINT, NSET=ALL\n";
    inp << "U\n";
    inp << "*EL PRINT, ELSET=LAYER1, FREQUENCY=0\n";
    inp << "  S\n";
    inp << "*END STEP\n";

    inp.close();

    // -------------------------------
    // Run CalculiX
    // -------------------------------
    std::string cmd = "ccx " + BASE_FILE_NAME + " > /dev/null 2>&1";
    int result = std::system(cmd.c_str());
    if (result != 0)
        throw std::runtime_error("CalculiX run failed");

    std::cout << "[FEA] CalculiX simulation complete with pressure = " << pressure << " MPa\n";


    std::ifstream output("calculix.dat");
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

    C3t3 deformed_mesh = mesh.get_c3t3();

    //------------------------------------------------------
    // Build id -> vertex-handle map *for the copy*
    // (iterate in the *same* order you used when writing *NODE)
    //------------------------------------------------------
    std::unordered_map<std::size_t, Tr::Vertex_handle> id_to_vh_copy;
    id = 1;
    for(auto vit = deformed_mesh.triangulation().finite_vertices_begin();
        vit != deformed_mesh.triangulation().finite_vertices_end(); ++vit, ++id)
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

std::cout << "[FEA] Applied " << id_to_displacement.size() 
          << " displacements to deformed mesh\n";

    return deformed_mesh;
}


void FEA::get_cut_stress_field(double pressure)
{
    const auto& nodes = mesh.get_nodes();
    const auto& c3t3 = mesh.get_c3t3();
    const auto& tr = c3t3.triangulation();

    const std::string BASE_FILE_NAME = "calculix";

    std::ofstream inp(BASE_FILE_NAME + ".inp");
    if (!inp.is_open())
        throw std::runtime_error("Failed to open CalculiX input file");

    // -------------------------------
    // *NODE block
    // -------------------------------
    inp << "*NODE\n";
    std::map<Tr::Vertex_handle, size_t> node_id_map;
    std::map<size_t, Tr::Vertex_handle> id_node_map;
    size_t id = 1;

    for (auto vit = tr.finite_vertices_begin(); vit != tr.finite_vertices_end(); ++vit, ++id) 
    {
        const Point p = vit->point().point();
        node_id_map[vit] = id;
        id_node_map[id] = vit;
        inp << id << ", " << p.x() << ", " << p.y() << ", " << p.z() << "\n";
    }

    // -------------------------------
    // *ELEMENT block (C3D10)
    // -------------------------------
    inp << "*ELEMENT, TYPE=C3D10, ELSET=TRENCH_ELEMS\n";
    size_t elem_id = 1;
    for (auto cit = tr.finite_cells_begin(); cit != tr.finite_cells_end(); ++cit, ++elem_id) 
    {
        inp << elem_id;
        for (int i = 0; i < 4; ++i)
            inp << ", " << node_id_map.at(cit->vertex(i));
        inp << "\n";
    }

    // -------------------------------
    // MATERIAL and SECTION
    // -------------------------------
    inp << "*MATERIAL, NAME=ALUMINUM\n";
    inp << "*ELASTIC\n";
    inp << "69000, 0.33\n"; // MPa
    inp << "*SOLID SECTION, ELSET=TRENCH_ELEMS, MATERIAL=ALUMINUM\n";

    // -------------------------------
    // BOUNDARY CONDITIONS
    // -------------------------------
    inp << "*BOUNDARY\n";
    for (const auto& node : nodes) 
    {
        if (node->is_boundary) 
        {
            inp << node_id_map.at(node->handle) << ", 1, 3, 0.0\n";
        }
    }

    // -------------------------------
    // SURFACE SET for trench wall
    // -------------------------------
    inp << "*SURFACE, TYPE=ELEMENT, NAME=TRENCH_FACES\n";
    std::unordered_set<size_t> trench_element_ids;
    size_t elem_counter = 1;

    for (auto cit = tr.finite_cells_begin(); cit != tr.finite_cells_end(); ++cit, ++elem_counter) {
        bool has_force_node = false;
        for (int i = 0; i < 4; ++i) {
            auto vh = cit->vertex(i);
            const auto& node = nodes[node_id_map.at(vh) - 1];
            if (node->status == NodeStatus::FORCE_NODE)
                has_force_node = true;
        }
        if (has_force_node) {
            trench_element_ids.insert(elem_counter);
            inp << elem_counter << ", S1\n"; // or S2/S3/S4 if you know which face is the wall
        }
    }

    // -------------------------------
    // LOAD STEP (apply pressure + output)
    // -------------------------------
    inp << "*STEP\n";
    inp << "*STATIC\n";
    inp << "*DLOAD\n";
    inp << "TRENCH_FACES, P, " << pressure << "\n";
    inp << "*NODE PRINT, NSET=ALL\nU\n";
    inp << "*EL FILE\nS\n";
    inp << "*END STEP\n";

    inp.close();

    // -------------------------------
    // Run CalculiX
    // -------------------------------
    std::string cmd = "ccx " + BASE_FILE_NAME + " > /dev/null 2>&1";
    int result = std::system(cmd.c_str());
    if (result != 0)
        throw std::runtime_error("CalculiX run failed");

    std::cout << "[FEA] CalculiX simulation complete with pressure = "
              << pressure << " MPa\n";

    // -------------------------------
    // Parse .frd for stress tensors
    // -------------------------------
    std::ifstream frd("calculix.frd");
    if (!frd.is_open())
        throw std::runtime_error("Cannot open calculix.frd");

    std::unordered_map<Tr::Vertex_handle, StressTensor> handle_to_tensor;
    std::string line;
    bool in_frame = false;

    while (std::getline(frd, line)) 
    {
        if (line.rfind("3C", 0) == 0) 
        {
            in_frame = !in_frame;
            continue;
        }
        if (!in_frame) continue;
        if (line.size() < 2 || line[0] != ' ' || line[1] != '1') continue;

        // Expected layout: " 1 <nodeID> 1  S  <Sxx> <Syy> <Szz> <Sxy> <Sxz> <Syz>"
        int rec, node_id, one;
        char key;
        double sxx, syy, szz, sxy, sxz, syz;
        std::istringstream iss(line);
        iss >> rec >> node_id >> one >> key >> sxx >> syy >> szz >> sxy >> sxz >> syz;

        auto it = id_node_map.find(node_id);
        if (it != id_node_map.end())
            handle_to_tensor[it->second] = {sxx, syy, szz};
    }

    // -------------------------------
    // Write back to mesh
    // -------------------------------
    for (auto& node : mesh.get_nodes()) 
    {
        auto it = handle_to_tensor.find(node->handle);
        if (it != handle_to_tensor.end())
            node->stress = it->second;
    }

    // -------------------------------
    // Export CSV
    // -------------------------------
    std::ofstream out("stress_field.csv");
    out << "x,y,z,sxx,syy,szz\n";
    for (const auto& node : mesh.get_nodes()) 
    {
        out << node->position.x() << ","
            << node->position.y() << ","
            << node->position.z() << ","
            << node->stress.xx << ","
            << node->stress.yy << ","
            << node->stress.zz << "\n";
    }
}


/*
void FEA::get_full_stress_field()
{

}
*/
