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


    std::map<Tr::Vertex_handle, size_t> node_id_map;
    std::unordered_set<size_t> corner_node_ids;
    std::map<std::pair<size_t, size_t>, size_t> midside_node_map;
    size_t id = 1;

    auto make_edge_key = [&](size_t a, size_t b) 
    {
        return (a < b) ? std::pair{a,b} : std::pair{b,a};
    };

    auto get_midside_node = [&](Tr::Vertex_handle a, Tr::Vertex_handle b) -> size_t
    {
        // get corner IDs
        size_t idA = node_id_map[a];
        size_t idB = node_id_map[b];

        std::pair<size_t, size_t> key = make_edge_key(idA, idB);

        // Already created?
        auto it = midside_node_map.find(key);
        if (it != midside_node_map.end())
            return it->second;

        // Create midpoint
        const auto& p1 = a->point().point();
        const auto& p2 = b->point().point();
        Point mid( (p1.x()+p2.x())*0.5,
                   (p1.y()+p2.y())*0.5,
                   (p1.z()+p2.z())*0.5 );

        // New node ID
        size_t new_id = id++;        // IMPORTANT: uses same node-id counter


        // Export midside node to .inp file
        inp << new_id << ", "
            << mid.x() << ", " << mid.y() << ", " << mid.z() << "\n";

        // Store midside ID
        midside_node_map[key] = new_id;
        return new_id;
    };

    // -------------------------------
    // *NODE block
    // -------------------------------
    inp << "*NODE\n";

    for (auto vit = tr.finite_vertices_begin(); 
         vit != tr.finite_vertices_end(); 
         ++vit)
    {
        const Point p = vit->point().point();
        node_id_map[vit] = id;
        corner_node_ids.insert(id);

        inp << id << ", " << p.x() << ", " << p.y() << ", " << p.z() << "\n";

        id++;   // increment ONLY here
    }

    for (auto cit = tr.finite_cells_begin(); 
         cit != tr.finite_cells_end(); 
         ++cit)
    {
        Tr::Vertex_handle v1 = cit->vertex(0);
        Tr::Vertex_handle v2 = cit->vertex(1);
        Tr::Vertex_handle v3 = cit->vertex(2);
        Tr::Vertex_handle v4 = cit->vertex(3);

        // Skip bad-quality tets just like you do for element writing
        if (tet_quality(v1->point().point(),
                        v2->point().point(),
                        v3->point().point(),
                        v4->point().point()) < 0.05)
            continue;

        // Generate midside nodes for all 6 edges
        get_midside_node(v1, v2);  // n12
        get_midside_node(v1, v3);  // n13
        get_midside_node(v1, v4);  // n14
        get_midside_node(v2, v3);  // n23
        get_midside_node(v2, v4);  // n24
        get_midside_node(v3, v4);  // n34
    }

    inp << "*NSET, NSET=ALL, GENERATE\n1, " << id-1 << ", 1\n";


    // -------------------------------
    // *ELEMENT block (C3D10)
    // -------------------------------
    //inp << "*ELEMENT, TYPE=C3D10, ELSET=TRENCH_ELEMS\n";
    //will switch to C3D10 later because it's better but C3D4 is easier to get this working first
    inp << "*ELEMENT, TYPE=C3D10, ELSET=TRENCH_ELEMS\n";

    std::unordered_map<Tr::Cell_handle, size_t> cell_to_elem_id;
    size_t elem_id = 1;
    for (auto cit = tr.finite_cells_begin(); 
         cit != tr.finite_cells_end(); 
         ++cit)
    {
        Tr::Vertex_handle v1 = cit->vertex(0);
        Tr::Vertex_handle v2 = cit->vertex(1);
        Tr::Vertex_handle v3 = cit->vertex(2);
        Tr::Vertex_handle v4 = cit->vertex(3);

        if (tet_quality(v1->point().point(),
                        v2->point().point(),
                        v3->point().point(),
                        v4->point().point()) < 0.05)
            continue;

        size_t n1 = node_id_map[v1];
        size_t n2 = node_id_map[v2];
        size_t n3 = node_id_map[v3];
        size_t n4 = node_id_map[v4];

        size_t n5  = midside_node_map[make_edge_key(n1,n2)]; // 1-2
        size_t n6  = midside_node_map[make_edge_key(n2,n3)]; // 2-3
        size_t n7  = midside_node_map[make_edge_key(n3,n1)]; // 3-1
        size_t n8  = midside_node_map[make_edge_key(n1,n4)]; // 1-4
        size_t n9  = midside_node_map[make_edge_key(n2,n4)]; // 2-4
        size_t n10 = midside_node_map[make_edge_key(n3,n4)]; // 3-4
        
        cell_to_elem_id[cit] = elem_id;

        inp << elem_id++ << ", "
            << n1  << ", " << n2  << ", " << n3  << ", " << n4  << ", "
            << n5  << ", " << n6  << ", " << n7  << ", "
            << n8  << ", " << n9  << ", " << n10 << "\n";
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

    for (auto cit = tr.finite_cells_begin(); cit != tr.finite_cells_end(); ++cit)
    {
        auto it_elem = cell_to_elem_id.find(cit);
        if (it_elem == cell_to_elem_id.end())
        {
            continue; // this cell was skipped (bad quality) → no element
        }

        size_t eid = it_elem->second;
        
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

                trench_faces.push_back({eid, face + 1});
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
                    in_displacement = true; continue;
            }

            if (in_displacement && line.empty())
            {
                    continue;
            }

            if (in_displacement)
            {
                    std::size_t id; double dx,dy,dz;

                    // Since CGAL can only represent C3D4 meshes we can only export corner nodes
                    /*
                    if (corner_node_ids.find(id) == corner_node_ids.end())
                    {
                        continue;
                    }
                    */

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
    const auto& tr   = c3t3.triangulation();

    std::string BASE_FILE_NAME = "calculix";

    std::ofstream inp(BASE_FILE_NAME + ".inp");
    if (!inp.is_open())
        throw std::runtime_error("Failed to open CalculiX input file");

    // -----------------------------------------------------
    // NODES
    // -----------------------------------------------------
    inp << "*NODE\n";
    std::map<Tr::Vertex_handle, size_t> node_id_map;
    size_t id = 1;

    for (auto vit = tr.finite_vertices_begin(); vit != tr.finite_vertices_end(); ++vit, ++id)
    {
        const Point p = vit->point().point();
        node_id_map[vit] = id;
        inp << id << ", " << p.x() << ", " << p.y() << ", " << p.z() << "\n";
    }

    inp << "*NSET, NSET=ALL, GENERATE\n1," << id-1 << ",1\n";

    // -----------------------------------------------------
    // ELEMENTS (C3D4 for now)
    // -----------------------------------------------------
    inp << "*ELEMENT, TYPE=C3D4, ELSET=TRENCH_ELEMS\n";
    size_t elem_id = 1;
    for (auto cit = tr.finite_cells_begin(); cit != tr.finite_cells_end(); ++cit, ++elem_id)
    {
        if (tet_quality(cit->vertex(0)->point().point(),
                        cit->vertex(1)->point().point(),
                        cit->vertex(2)->point().point(),
                        cit->vertex(3)->point().point()) < 0.05)
            continue;

        inp << elem_id;
        for (int i = 0; i < 4; ++i)
            inp << ", " << node_id_map.at(cit->vertex(i));
        inp << "\n";
    }

    // -----------------------------------------------------
    // MATERIAL
    // -----------------------------------------------------
    inp << "*MATERIAL, NAME=ALUMINUM\n";
    inp << "*ELASTIC\n69000., 0.33\n";
    inp << "*SOLID SECTION, ELSET=TRENCH_ELEMS, MATERIAL=ALUMINUM\n";

    // -----------------------------------------------------
    // BOUNDARIES
    // -----------------------------------------------------
    inp << "*BOUNDARY\n";
    for (const auto& node : nodes)
        if (node->is_boundary)
            inp << node_id_map.at(node->handle) << ", 1, 3, 0.0\n";

    // -----------------------------------------------------
    // SURFACE & PRESSURE (same trench face logic)
    // -----------------------------------------------------
    std::vector<FaceRef> trench_faces;
    size_t elem_counter = 1;

    for (auto cit = tr.finite_cells_begin(); cit != tr.finite_cells_end(); ++cit, ++elem_counter)
    {
        for (int face = 0; face < 4; ++face)
        {
            auto vh0 = cit->vertex((face + 1) % 4);
            auto vh1 = cit->vertex((face + 2) % 4);
            auto vh2 = cit->vertex((face + 3) % 4);

            const auto& n0 = nodes[node_id_map.at(vh0) - 1];
            const auto& n1 = nodes[node_id_map.at(vh1) - 1];
            const auto& n2 = nodes[node_id_map.at(vh2) - 1];

            if (n0->status == NodeStatus::FORCE_NODE ||
                n1->status == NodeStatus::FORCE_NODE ||
                n2->status == NodeStatus::FORCE_NODE)
            {
                trench_faces.push_back({elem_counter, face + 1});
            }
        }
    }

    inp << "*SURFACE, TYPE=ELEMENT, NAME=TRENCH_FACES\n";
    for (auto& f : trench_faces)
        inp << f.elem_id << ", S" << f.face_id << "\n";

    // -----------------------------------------------------
    // LOAD STEP
    // -----------------------------------------------------
    inp << "*STEP\n*STATIC\n*DLOAD\n";
    inp << "TRENCH_FACES, P, " << pressure << "\n";
    inp << "*NODE FILE, NSET=ALL\nU\n";
    inp << "*EL FILE, ELSET=TRENCH_ELEMS\nS\n";
    inp << "*NODE PRINT, NSET=ALL\nU\n";
    inp << "*EL PRINT, ELSET=TRENCH_ELEMS\nS\n";
    inp << "*END STEP\n";
    inp.close();

    // -----------------------------------------------------
    // Run CalculiX
    // -----------------------------------------------------
    if (std::system(("ccx " + BASE_FILE_NAME + " > /dev/null 2>&1").c_str()) != 0)
        throw std::runtime_error("CalculiX run failed");

    std::cout << "[FEA] Simulation complete, extracting stress field.\n";

    // -----------------------------------------------------
    // Parse calculix.dat for stress tensors
    // -----------------------------------------------------
    std::ifstream dat("calculix.dat");
    if (!dat.is_open())
        throw std::runtime_error("Cannot open calculix.dat");

    // -----------------------------------------------------
    // 1. Parse stresses from CalculiX .dat file
    // -----------------------------------------------------
    std::unordered_map<size_t, StressTensor> id_to_stress;
    std::string line;
    bool in_stress = false;

    while (std::getline(dat, line))
    {
        // Find the start of the stress block
        if (line.find("stresses") != std::string::npos)
        {
            std::cout << "found stresses" << std::endl;
            in_stress = true;
            continue;
        }

        // End of block (empty line or another header)
        if (in_stress && (line.empty() || line[0] == '-'))
        {
            continue;
        }

        if (!in_stress) continue;

        std::istringstream iss(line);
        size_t elem_id;
        int ip; // integration point (ignore for C3D4)
        double sxx, syy, szz, sxy, sxz, syz;

        // Parse 8 fields: element, integ.pnt., 6 stresses
        if (iss >> elem_id >> ip >> sxx >> syy >> szz >> sxy >> sxz >> syz)
        {
            // Store only the normal components for now
            id_to_stress[elem_id] = StressTensor{sxx, syy, szz};
        }
    }

    dat.close();

    std::cout << "[FEA] Parsed " << id_to_stress.size() << " element stresses from .dat\n";

    // -----------------------------------------------------
    // 2. Average stresses to nodes
    // -----------------------------------------------------
    std::unordered_map<Tr::Vertex_handle, StressTensor> avg_stress;
    std::unordered_map<Tr::Vertex_handle, int> counts;

    // If you generated the .inp sequentially, element IDs = iteration order
    size_t eid = 1;
    for (auto cit = tr.finite_cells_begin(); cit != tr.finite_cells_end(); ++cit, ++eid)
    {

        if (tet_quality(cit->vertex(0)->point().point(),
                        cit->vertex(1)->point().point(),
                        cit->vertex(2)->point().point(),
                        cit->vertex(3)->point().point()) < 0.05)
            continue;

        auto it = id_to_stress.find(eid);
        if (it == id_to_stress.end()) continue;

        // Distribute the element stress to each vertex of the tetrahedron
        for (int i = 0; i < 4; ++i)
        {
            auto vh = cit->vertex(i);
            avg_stress[vh].xx += it->second.xx;
            avg_stress[vh].yy += it->second.yy;
            avg_stress[vh].zz += it->second.zz;
            counts[vh]++;
        }
    }

    // -----------------------------------------------------
    // 3. Export averaged nodal stresses to CSV
    // -----------------------------------------------------
    std::ofstream out("stress_field.csv");
    out << "x,y,z,sxx,syy,szz\n";

    for (auto node : mesh.get_nodes())
    {
        auto vh = node->handle;
        auto it = counts.find(vh);
        if (it == counts.end()) 
        {
            //std::cout << "Couldn't find vertex handle in counts" << std::endl;
            continue;
        }

        if (node->status == NodeStatus::CUT_NEXT || node->status == NodeStatus::FORCE_NODE)
        {
            double n = static_cast<double>(it->second);
            const auto& pt = vh->point();
            const auto& s = avg_stress[vh];

            out << pt.x() << "," << pt.y() << "," << pt.z() << ","
                << s.xx / n << "," << s.yy / n << "," << s.zz / n << "\n";
        }

    }

    out.close();
    std::cout << "[FEA] Stress field written to stress_field.csv\n";
}


/*
void FEA::get_full_stress_field()
{

}
*/
